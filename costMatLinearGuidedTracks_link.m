function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
    errFlag] = costMatLinearGuidedTracks_link(movieInfo,kalmanFilterInfoFrame1,...
    costMatParam,nnDistFeatures,probDim,prevCost,featLifetime,...
    trackedFeatureIndx,currentFrame)
%
% costMatLinearGuidedTracks_link: cost matrix for guided ('on rails') linear motion with stopping.
% To be used with u-track (http://lccb.hms.harvard.edu/software.html). This code is heavily
% influenced by plusTipCostMatLinearMotionLink.m by Khuloud Jaqaman.
% 
% Fabian Peters 2015 
%
%INPUT
%       movieInfo             : An nx1 array (n = number of frames in
%                               movie) containing the fields:
%             .allCoord           : x,dx,y,dy,[z,dz] of features collected in one
%                                   matrix.
%             .amp                : Amplitudes of PSFs fitting detected features.
%                                   1st column for values and 2nd column
%                                   for standard deviations.
%             .num                : Number of features in each frame.
%             .nnDist             : Distance from each feature to its nearest
%                                   neighbor. Not needed at the moment.
%
%      kalmanFilterInfoFrame1 : Structure with at least the following fields:
%             .stateVec           : Kalman filter state vector for each
%                                   feature in 1st frame.
%             .stateCov           : Kalman filter state covariance matrix
%                                   for each feature in 1st frame.
%             .noiseVar           : Variance of state noise for each
%                                   feature in 1st frame.
%
%      costMatParam           : Structure with fields:
%             .linearMotion       : 1 to propagate enforcing a linear motion
%                                   model (no stopping), 0 otherwise.
%             .minSearchRadius    : Minimum allowed search radius for brownStdMult search window
%             .maxSearchRadius    : Maximum allowed search radius for brownStdMult search window
%             .brownStdMult       : Somewhat unfortunately named. Factor multiplied with noise
%                                   (std) from kalman filter to get search radius.
%                                   Overruled by maxSpeed.
%             .useLocalDensity    : Logical variable indicating whether to use
%                                   local density in brownStdMult search radius estimation.
%             .nnWindow           : Number of past frames for calculating
%                                   nearest neighbor distance in brownStdMult search radius estimation.
%             .maxSpeed           : Maximum displacement between two frames.
%             .maxVelocityAngle   : Max angle between current velocity (as indicated by kalman
%                                   filter) and vector connecting current and new position. Only
%                                   applies if dist > minSpeedAngleFilter.
%             .maxYdistSlow       : Max Y displacement between two frames if dist < minSpeedAngleFilter.
%             .maxYdistFast       : Max Y displacement between two frames if dist >= minSpeedAngleFilter.
%             .minSpeedAngleFilter: Min distance in px from which maxVelocityAngle and
%                                   maxHorizontalAngle are applied. Helps to ignore small drifts. 
%             .maxAmpRatio        : Max acceptable amp ratio for the linking of two particles.
%             .distFact,ampFact   : weight of distance and amp difference in calculation of total
%                                   cost. Total cost is calculated as
%                                   distFact*(distance/maxSpeed)^2 + ampFact*(ampCost/maxAmpRatio)
%
%      nnDistFeatures         : Matrix of nearest neighbor distances of
%                               features in first frame as well as of
%                               features in previous frames that they are
%                               connected to. NOT USED.
%      probDim                : Problem dimensionality. 2 (for 2D) or 3 (for 3D). NEEDS TO BE 2.
%      prevCost               : Matrix of previous linking costs.
%      featLifetime           : Lengths of tracks that features in
%                               first frame belong to.  NOT USED.
%      trackedFeatureIndx     : The matrix of feature index connectivity up
%                               to current frame.  NOT USED.
%      currentFrame           : Current frame that is being linked to the
%                               next frame.
%
%OUTPUT
%       costMat               : Cost matrix.
%       propagationScheme     : Propagation scheme corresponding to each
%                               cost in the cost matrix.
%       kalmanFilterInfoFrame2: Structure with at least the following fields:
%             .stateVec           : Kalman filter prediction of state
%                                   vector in 2nd frame based on all 3
%                                   motion models.
%             .stateCov           : Kalman filter prediction of state
%                                   covariance matrix in 2nd frame based on
%                                   all 3 motion models.
%             .obsVec             : Kalman filter prediction of the
%                                   observed variables in 2nd frame based
%                                   on all 3 motion models.
%       nonlinkMarker         : Value indicating that a link is not allowed.
%       errFlag               : 0 if function executes normally, 1 otherwise.
%
%
% you can redistribute and/or modify this file under the terms of the GNU General Public License
% as published by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 

%% initial checks

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatLinearGuidedTracks_link')
    disp('--costMatLinearGuidedTracks_link: Incorrect number of input arguments.');
    errFlag  = 1;
    return
end

% ensure that the problem is two-dimensional (not tested for other problems)
if probDim ~= 2
    disp('--costMatLinearGuidedTracks_link: Problem dimension must be 2.');
    errFlag  = 1;
    return
end


%% Output

costMat = [];
propagationScheme = [];
kalmanFilterInfoFrame2 = [];
nonlinkMarker = [];
errFlag = 0;

%% Input

linearMotion = costMatParam.linearMotion;
maxSpeed = costMatParam.maxSpeed;
maxYdistFast = costMatParam.maxYdistFast;
maxYdistSlow = costMatParam.maxYdistSlow;
brownStdMult = costMatParam.brownStdMult;
useLocalDensity = costMatParam.useLocalDensity;
nnWindow = costMatParam.nnWindow;
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
maxVelocityAngle = costMatParam.maxVelocityAngle;
minSpeedAngleFilter = costMatParam.minSpeedAngleFilter;
maxAmpRatio = costMatParam.maxAmpRatio;
distFact = costMatParam.distFact;
ampFact = costMatParam.ampFact;

%extract the two frames of interest from movieInfo
movieInfo = movieInfo(currentFrame:currentFrame+1);

%% Motion propagation

% number of propagation schemes used:
% one for forward movement and one for stopping
numSchemes = 2;

%calculate vector sizes
vecSize = 2 * probDim;

%construct transition matrices
if linearMotion
    transMat(:,:,1) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
else
    transMat(:,:,1) = eye(vecSize); %zero drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix    
end

%construct observation matrix
observationMat = [eye(probDim) zeros(probDim)]; %observation matrix

%get number of features in the 2 frames
numFeaturesFrame1 = movieInfo(1).num;
numFeaturesFrame2 = movieInfo(2).num;

%reserve memory for "kalmanFilterInfoframe2"
kalmanFilterInfoFrame2 = struct('stateVec',zeros(numFeaturesFrame1,vecSize,numSchemes),...
    'stateCov',zeros(vecSize,vecSize,numFeaturesFrame1,numSchemes),...
    'obsVec',zeros(numFeaturesFrame1,probDim,numSchemes));

% store previous velocity for each feature to disable linking at wrong angles
oldVelocity = zeros(numFeaturesFrame1,2);
% store the actual last position of each feature
oldCoord = zeros(numFeaturesFrame1,2);

%apply Kalman filters to each feature in 1st frame
for iFeature = 1 : numFeaturesFrame1
    
    %get state vector and its covariance matrix of feature in 1st frame
    stateOld = kalmanFilterInfoFrame1.stateVec(iFeature,:)';
    stateCovOld = kalmanFilterInfoFrame1.stateCov(:,:,iFeature);
    noiseVar = abs(kalmanFilterInfoFrame1.noiseVar(:,:,iFeature));
    
    oldVelocity(iFeature,:) = stateOld(end-1:end);
    oldCoord(iFeature,:) = stateOld(1:2);
    
    %go over all possible propagation schemes
    for iScheme = 1 : numSchemes
        
        %predict state vector of feature in 2nd frame
        stateVec = transMat(:,:,iScheme)*stateOld;
        
        %predict state covariance matrix of feature in 2nd frame
        stateCov = transMat(:,:,iScheme)*stateCovOld*transMat(:,:,iScheme)' ...
            + noiseVar;
        
        %determine observation vector of feature in 2nd frame (i.e. the
        %propagated position of the feature)
        obsVec = observationMat*stateVec;
        
        %save information in kalmanFilterInfoFrame2
        kalmanFilterInfoFrame2.stateVec(iFeature,:,iScheme) = stateVec';
        kalmanFilterInfoFrame2.stateCov(:,:,iFeature,iScheme) = stateCov;
        kalmanFilterInfoFrame2.obsVec(iFeature,:,iScheme) = obsVec';
        
    end
    
end

%get the propagated positions of features in 1st frame based on the two propagation schemes
propagatedPos = kalmanFilterInfoFrame2.obsVec;

%put the coordinates of features in the 2nd frame in one matrix
coord2 = movieInfo(2).allCoord(:,1:2:end);

distCostMat = zeros(numFeaturesFrame1, numFeaturesFrame2, numSchemes);

%calculate the cost matrices for all propagation schemes
for iScheme = 1 : numSchemes
    
    %put the propagated x and y coordinates of features from 1st frame in
    %one matrix
    coord1 = propagatedPos(:,:,iScheme);
    
    %calculate the distances between features and their projected positions
    distCostMat(:,:,iScheme) = createDistanceMatrix(coord1,coord2);    
end

%find the minimum cost for the link between every pair, which also
%determines the best propagation scheme to perform that link
[distCostMat,propagationScheme] = min(distCostMat,[],3);

%% apply distance limits / max speed

% calculate distances between features and their true previous positions
trueDistMat = createDistanceMatrix(oldCoord,coord2);
distCostMat(trueDistMat > maxSpeed) = NaN;

% calculate distance in y direction only
yDistMat = createDistanceMatrix([zeros(size(oldCoord,1),1) oldCoord(:,2)], ...
    [zeros(size(coord2,1),1) coord2(:,2)]);

% limit possible links according to max y distance
distCostMat((distCostMat < minSpeedAngleFilter) & (yDistMat > maxYdistSlow)) = NaN;
distCostMat((distCostMat >= minSpeedAngleFilter) & (yDistMat > maxYdistFast)) = NaN;

%% Limit search radius according to previous velocity

if useLocalDensity
    closestDistScale = 2;
    maxStdMult = 100;
end

%calculate nearest neighbor distance given feature history
frameNum = size(nnDistFeatures,2);
tmpNN = max(1,frameNum-nnWindow);
nnDistTracks = min(nnDistFeatures(:,tmpNN:end),[],2);

%determine which features are not first appearances
notFirstAppearance = squeeze(kalmanFilterInfoFrame1.noiseVar(1,1,:)) >= 0;

%get the Kalman standard deviation of all features in frame 1
kalmanStd = sqrt(probDim * abs(squeeze(kalmanFilterInfoFrame1.noiseVar(1,1,:))));

%copy brownStdMult into vector
stdMultInd = repmat(brownStdMult,numFeaturesFrame1,1);

%if local density information is used to expand search radius ...
if useLocalDensity

    %divide each feature's nearest neighbor distance/closestDistScale by kalmanStd
    ratioDist2Std = nnDistTracks./kalmanStd/closestDistScale;

    %make ratios larger than maxStdMult equal to maxStdMult
    ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

    %expand search radius multiplication factor if possible
    stdMultInd = max([stdMultInd ratioDist2Std],[],2);

end

%get the search radius of each feature in frame 1 and make sure it falls
%within reasonable limits
searchRadius = stdMultInd .* kalmanStd;
searchRadius((searchRadius>maxSearchRadius)&notFirstAppearance) = maxSearchRadius;
searchRadius((searchRadius<minSearchRadius)&notFirstAppearance) = minSearchRadius;

%replicate the search radius to compare to cost matrix
searchRadius = repmat(searchRadius,1,numFeaturesFrame2);

%assign NaN to costs corresponding to distance > searchRadius
distCostMat(trueDistMat > searchRadius) = NaN;

%% apply velocity angle limits
% these are only applied if the particle is moving at a minimum speed in order to avoid breaking
% tracks due to slight shifts of the particle centre

%calculate angles between the two velocities
velocityAngleMat = zeros(numFeaturesFrame1, numFeaturesFrame2);
horizontalAngleMat = zeros(numFeaturesFrame1, numFeaturesFrame2);
for iF1=1:numFeaturesFrame1
    oldVel = oldVelocity(iF1,:);
    for iF2=1:numFeaturesFrame2
        % the computations below are expensive; only perform them if the particle is still a
        % candidate
        if isnan(distCostMat(iF1,iF2))
           velocityAngleMat(iF1,iF2) = NaN;
           horizontalAngleMat(iF1,iF2) = NaN;
           continue;
        end
        
        % velocity of the connection tested
        newVel = coord2(iF2,:)-oldCoord(iF1,:);
        
        % some angles are NaN because either velocity is a zero vector. Set
        % these to 0
        if norm(newVel)==0 || norm(oldVel)==0
            velocityAngle = 0;
        else
            velocityAngle = vvAngle(oldVel, newVel);    
        end        
        velocityAngleMat(iF1, iF2) = velocityAngle;
        
        % angle between new velocity and the horizontal line
        horizontalAngle = vvAngle([1 0], newVel);
        if horizontalAngle>90,
            horizontalAngle = 180 - horizontalAngle;
        end        
        % a lot of angles are NaN because the old velocity is a zero vector. Set
        % these to 0
        if isnan(horizontalAngle)
            horizontalAngle = 0;
        end
        horizontalAngleMat(iF1, iF2) = horizontalAngle;
    end
end

% treat angle limits differently depending on whether the particle in frame 1 is  a new appearance
% or already part of a track
notFirstAppearanceMatrix = repmat(notFirstAppearance, [1 numFeaturesFrame2]);

% ignore angles for particles which were previously slower than the angle requirement
oldVelocityMag = zeros(size(oldVelocity,1),1);
for iV=1:numel(oldVelocityMag)
    oldVelocityMag(iV) = norm(oldVelocity(iV,:));
end
oldVelocityMagMatrix = repmat(oldVelocityMag, [1 numFeaturesFrame2]);

% limit search / distCostMatrix according to maxVelocityAngle, given that the particle from frame 1 
% has sufficient velocity history
distCostMat((velocityAngleMat > maxVelocityAngle) &...
    (distCostMat >= minSpeedAngleFilter) &...
    (notFirstAppearanceMatrix == 1) &...
    (oldVelocityMagMatrix >= minSpeedAngleFilter)) = NaN;

% limit new particles to maxVelocityAngle relative to the horizontal line
distCostMat((horizontalAngleMat > maxVelocityAngle) &...
    (distCostMat >= minSpeedAngleFilter) &...
    ((notFirstAppearanceMatrix == 0) |...
    (oldVelocityMagMatrix < minSpeedAngleFilter))) = NaN;


%% Amplitude factor

%put feature amplitudes from both frames in vectors
amp1 = movieInfo(1).amp(:,1);
amp2 = movieInfo(2).amp(:,1);

%make a matrix of amplitude costs
ampRatioMat = repmat(amp1,1,numFeaturesFrame2)./repmat(amp2',numFeaturesFrame1,1);
ampRatioMat(ampRatioMat<1) = 1./ampRatioMat(ampRatioMat<1);

% exclude overly big differences
ampRatioMat(ampRatioMat > maxAmpRatio) = NaN;


%% combine costs

%square the cost matrix to make the cost = distance squared
distCostMat = distCostMat.^2;
% normalize
distCostMat = distCostMat/maxSpeed^2;
ampCost = ampRatioMat/maxAmpRatio;

% add using factors
costMat = distFact*distCostMat + ampFact*ampCost;


%% Birth and death

%append matrix to allow birth and death
if isstruct(prevCost)
    prevCostMax = prevCost.max;
else
    prevCostMax = max(prevCost(:));
end

if ~isnan(prevCostMax) && prevCostMax ~= 0
    maxCost = 1.05*prevCostMax;
else
    maxCost = max(prctile(costMat(:),80),eps);
end
deathCost = maxCost * ones(numFeaturesFrame1,1);
birthCost = maxCost * ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);
lrBlock = costMat';
lrBlock(~isnan(lrBlock)) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];


%% nonLinkMarker

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;   
    
end

