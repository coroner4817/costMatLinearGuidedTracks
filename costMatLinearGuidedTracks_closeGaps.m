function [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,errFlag] = ...
        costMatLinearGuidedTracks_closeGaps(trackedFeatInfo,trackedFeatIndx,trackStartTime,...
        trackEndTime,costMatParam,gapCloseParam,kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)

%
% To be used with u-track (http://lccb.hms.harvard.edu/software.html).
% costMatLinearGuidedTracks_closeGaps provides a cost matrix for closing gaps and capturing merges/splits
% in tracks linked by costMatLinearGuidedTracks_link. It is based on costMatRandomDirectedSwitchingMotionLink
% by Khuloud Jaqaman and plusTipCostMatCloseGaps by Kathryn Applegate
% 
% Fabian Peters 2015 
%
%INPUT  trackedFeatInfo: The positions and amplitudes of the tracked
%                        features from linkFeaturesKalman.
%                        Number of rows = number of tracks.
%                        Number of columns = 8*number of frames.
%                        Each row consists of
%                        [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                        in image coordinate system (coordinates in
%                        pixels). NaN is used to indicate time points
%                        where the track does not exist.
%       trackedFeatIndx: Connectivity matrix of features between frames.
%                        Rows indicate continuous tracks, while columns
%                        indicate frames. A track that ends before the
%                        last time point is followed by zeros, and a track
%                        that starts at a time after the first time point
%                        is preceded by zeros.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       costMatParam   : Structure containing variables needed for cost
%                        calculation. Contains the fields:
%          .fluctRad             : Size in pixels of tube radius around track
%                                trajectory. The search region in the backward
%                                direction will expand out from the track at
%                                the final point.  This value also determines
%                                the search radius around the final point
%                                wherein any candidates will be considered for
%                                forward linking, even if they fall slightly
%                                behind the point.  This ensures that tracks
%                                starting from a fluctuation during a pause
%                                will still be picked up as candidates for
%                                pause.
%          .maxFAngle          : Max angle in degrees allowed between the end
%                                track's final velocity vector and the
%                                displacement vector between end and start.
%                                Also the max angle between the end and start
%                                tracks themselves.
%          .maxBAngle          : Angle in degrees used to expand backward
%                                search region, giving a distance-dependent
%                                criterion for how far a start track
%                                can be from the lattice to be considered a
%                                candidate for linking. THIS IS CURRENTLY A
%                                HARDWIRED PARAMETER
%          .backVelMultFactor : Muliplication factor of max growth speed used
%                               to define candidate search area in the
%                               backward direction.
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to
%                           it.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           and 0 if merging and splitting are not allowed.
%                           For MT tracking, there are no merges/splits, so
%                           this should be 0.
%       kalmanFilterInfo: Structure array with number of entries equal to
%                         number of frames in movie. Contains the fields:
%             .stateVec   : Kalman filter state vector for each
%                           feature in frame.
%             .stateCov   : Kalman filter state covariance matrix
%                           for each feature in frame.
%             .noiseVar   : Variance of state noise for each
%                           feature in frame.
%             .stateNoise : Estimated state noise for each feature in
%                           frame.
%             .scheme     : 1st column: propagation scheme connecting
%                           feature to previous feature. 2nd column:
%                           propagation scheme connecting feature to
%                           next feature.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%       movieInfo      : movieInfo as input to trackCloseGapsKalman. Not
%                        really used in this code, but needed for
%                        compatibility with other cost functions.
%
%OUTPUT costMat       : Cost matrix.
%       nonlinkMarker : Value indicating that a link is not allowed.
%       indxMerge     : Index of tracks that have possibly merged with
%                       tracks that end before the last time points.
%       numMerge      : Number of such tracks.
%       indxSplit     : Index of tracks from which tracks that begin after
%                       the first time point might have split.
%       numSplit      : Number of such tracks.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%

%% Output

costMat = [];
nonlinkMarker = [];
indxMerge = [];
numMerge = [];
indxSplit = [];
numSplit = [];
errFlag = 0;


%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatLinearGuidedTracks_closeGaps')
    disp('--costMatLinearGuidedTracks_closeGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

% ensure that the problem is two-dimensional (not tested for other problems)
if probDim ~= 2
    disp('--costMatLinearGuidedTracks_link: Problem dimension must be 2.');
    errFlag  = 1;
    return
end

% merging / splitting currently not supported
if gapCloseParam.mergeSplit ~= 0
    disp('--costMatLinearGuidedTracks_link: merging / splitting currently not supported.');
    errFlag  = 1;
    return
end

%get cost matrix parameters
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
maxSpeed = costMatParam.maxSpeed;
brownStdMult = costMatParam.brownStdMult;
linStdMult   = costMatParam.linStdMult;
linScaling = costMatParam.linScaling;
timeReachConfL = costMatParam.timeReachConfL;
sin2AngleMax = (sin(costMatParam.maxAngleVV*pi/180))^2;
sin2AngleMaxVD = 0.5;
if isfield(costMatParam,'ampRatioLimit') && ~isempty(costMatParam.ampRatioLimit)
    minAmpRatio = costMatParam.ampRatioLimit(1);
    maxAmpRatio = costMatParam.ampRatioLimit(2);
    useAmp = 1;
else
    minAmpRatio = 0;
    maxAmpRatio = Inf;
    useAmp = 0;
end
if isfield(costMatParam,'gapPenalty') && ~isempty(costMatParam.gapPenalty)
    gapPenalty = costMatParam.gapPenalty;
else
    gapPenalty = 1;
end

%get gap closing parameters
timeWindow = gapCloseParam.timeWindow;

%make sure that timeReachConfB is <= timeWindow
timeReachConfL = min(timeReachConfL,timeWindow);

%find the number of tracks to be linked and the number of frames in the movie
[numTracks,numFrames] = size(trackedFeatInfo);
numFrames = numFrames / 8;

%list the tracks that start and end in each frame
tracksPerFrame = repmat(struct('starts',[],'ends',[]),numFrames,1);
for iFrame = 1 : numFrames    
    tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
    tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
end

%% Pre-processing

%get the x,y-coordinates and amplitudes at the starts of tracks
coordStart = zeros(numTracks,probDim);
ampStart   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordStart(iTrack,:) = full(trackedFeatInfo(iTrack,...
        (trackStartTime(iTrack)-1)*8+1:(trackStartTime(iTrack)-1)*8+probDim));
    ampStart(iTrack) = full(trackedFeatInfo(iTrack,(trackStartTime(iTrack)-1)*8+4));
end

%get the x,y-coordinates and amplitudes at the ends of tracks
coordEnd = zeros(numTracks,probDim);
ampEnd   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordEnd(iTrack,:) = full(trackedFeatInfo(iTrack,...
        (trackEndTime(iTrack)-1)*8+1:(trackEndTime(iTrack)-1)*8+probDim));
    ampEnd(iTrack) = full(trackedFeatInfo(iTrack,(trackEndTime(iTrack)-1)*8+4));
end

%% from estimTrackTypeParamRDS

%get number of tracks from initial linking and number of frames
[numTracksLink,numFrames] = size(trackedFeatIndx);

%reserve memory for output variables
xyzVelS = zeros(numTracksLink,probDim);
xyzVelE = zeros(numTracksLink,probDim);
noiseStd = zeros(numTracksLink,1);
trackMeanDispMag = NaN(numTracksLink,1);

%get the start times, end times and lifetimes of all tracks
trackSEL = getTrackSEL(trackedFeatInfo);
trackStartTime = trackSEL(:,1);
trackEndTime   = trackSEL(:,2);

%go over all tracks
for iTrack = 1 : numTracksLink
    
    %get current track's coordinates
    currentTrack = (reshape(trackedFeatInfo(iTrack,:)',8,[]))';
    currentTrack = currentTrack(:,1:probDim);
    currentTrack = full(currentTrack(trackStartTime(iTrack):trackEndTime(iTrack),:));
    
    %calculate the track's mean displacement
    if size(currentTrack,1) > 1
        trackMeanDispMag(iTrack) = mean(sqrt(sum(diff(currentTrack,1,1).^2,2)));
    end
    
    %assign velocity
    xyzVelS(iTrack,:) = kalmanFilterInfo(trackStartTime(...
        iTrack)).stateVec(trackedFeatIndx(iTrack,...
        trackStartTime(iTrack)),probDim+1:2*probDim);
    xyzVelE(iTrack,:) = kalmanFilterInfo(trackEndTime(...
        iTrack)).stateVec(trackedFeatIndx(iTrack,...
        trackEndTime(iTrack)),probDim+1:2*probDim);

    %assign noise std
    noiseStd(iTrack) = sqrt( abs( kalmanFilterInfo(trackEndTime(...
        iTrack)).noiseVar(1,1,trackedFeatIndx(iTrack,...
        trackEndTime(iTrack))) ) );           
end

%%

%calculate the average mean displacement for all tracks, to assign to
%tracks that have no mean displacement estimate
meanDispAllTracks = nanmean(trackMeanDispMag);


%% from getSearchRegionRDS

%reserve memory for output
longVecSAll  = zeros(probDim,timeWindow,numTracks);
longVecEAll  = zeros(probDim,timeWindow,numTracks);
longRedVecSAll  = zeros(probDim,timeWindow,numTracks);
longRedVecEAll  = zeros(probDim,timeWindow,numTracks);
shortVecSAll = zeros(probDim,timeWindow,numTracks);
shortVecEAll = zeros(probDim,timeWindow,numTracks);

%define square root of "problem dimension" to avoid calculating it many times
sqrtDim = sqrt(probDim);

%put time scaling of forward linear motion in a vector
timeScalingLin = [(1:timeReachConfL).^linScaling(1) ...
    (timeReachConfL)^linScaling(1) * (2:timeWindow-timeReachConfL+1).^linScaling(2)];
for iTrack = 1 : numTracks
          
    %get velocity, its magnitude and "direction of motion"
    %at track start
    velDriftS = xyzVelS(iTrack,:)';
    velMagS = sqrt(velDriftS' * velDriftS);
    directionMotionS = velDriftS / velMagS;
    %at track end
    velDriftE = xyzVelE(iTrack,:)';
    velMagE = sqrt(velDriftE' * velDriftE);
    directionMotionE = velDriftE / velMagE;

    %obtain vector(s) perpendicular to direction of motion
    %at track start
    perpendicularS = [-directionMotionS(2) directionMotionS(1)]';
    %at track end
    perpendicularE = [-directionMotionE(2) directionMotionE(1)]';

    %calculate the expected displacement due to drift for all time
    %gaps
    %at track start
    dispDrift1FS = velMagS * timeScalingLin;
    %at track end
    dispDrift1FE = velMagE * timeScalingLin;

    %calculate the expected displacement along x (= along y, [z]) due to
    %brownian motion for all time gaps
    dispBrown1 = brownStd(iTrack) * timeScalingBrown;

    %copy brownStdMult into vector that might be modified using
    %local density
    brownStdMultModS = brownStdMult'; %for track start
    brownStdMultModE = brownStdMult'; %for track end

    %determine the "long vectors" of the search rectangles for all time
    %gaps when direction of motion is continued
    %at track start
    longVec1FS = ...
        (repmat((linStdMult' .* dispDrift1FS),probDim,1) + ...
        repmat((brownStdMult' .* dispBrown1 * sqrtDim),probDim,1)) .* ...
        repmat(directionMotionS,1,timeWindow);
    longVecFSMag = sqrt((diag(longVec1FS' * longVec1FS))');  %magnitude
    longVecFSDir = longVec1FS ./ repmat(longVecFSMag,probDim,1); %direction
    %at track end
    longVec1FE = ...
        (repmat((linStdMult' .* dispDrift1FE),probDim,1) + ...
        repmat((brownStdMult' .* dispBrown1 * sqrtDim),probDim,1)) .* ...
        repmat(directionMotionE,1,timeWindow);
    longVecFEMag = sqrt((diag(longVec1FE' * longVec1FE))');  %magnitude
    longVecFEDir = longVec1FE ./ repmat(longVecFEMag,probDim,1); %direction

    %determine the "long vectors" of the search rectangles for all time
    %gaps when direction of motion is reversed
    %at start
    longVec1BS = ...
        repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
        repmat(directionMotionS,1,timeWindow);
    longVecBSMag = sqrt((diag(longVec1BS' * longVec1BS))');  %magnitude
    longVecBSDir = longVec1BS ./ repmat(longVecBSMag,probDim,1); %direction
    %at end
    longVec1BE = ...
        repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
        repmat(directionMotionE,1,timeWindow);
    longVecBEMag = sqrt((diag(longVec1BE' * longVec1BE))');  %magnitude
    longVecBEDir = longVec1BE ./ repmat(longVecBEMag,probDim,1); %direction

    %determine the "short vectors"
    %at track starts
    shortVecS1 = ...
        repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
        repmat(perpendicularS,1,timeWindow);
    shortVecSMag = sqrt((diag(shortVecS1' * shortVecS1))');  %magnitude
    shortVecSDir = shortVecS1 ./ repmat(shortVecSMag,probDim,1); %direction
    %at track ends
    shortVecE1 = ...
        repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
        repmat(perpendicularE,1,timeWindow);
    shortVecEMag = sqrt((diag(shortVecE1' * shortVecE1))');  %magnitude
    shortVecEDir = shortVecE1 ./ repmat(shortVecEMag,probDim,1); %direction

    %make sure that "long vectors" are longer than minimum allowed
    %at start
    longVecSMagTmp = max([longVecFSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
    longVec1FS = repmat(longVecSMagTmp,probDim,1) .* longVecFSDir; %new long vector
    %at end
    longVecEMagTmp = max([longVecFEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
    longVec1FE = repmat(longVecEMagTmp,probDim,1) .* longVecFEDir; %new long vector

    %make sure that backwards "long vectors" are
    %within allowed range
    %at start
    longVecSMagTmp = max([longVecBSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
    longVecSMagTmp = min([longVecSMagTmp;maxSearchRadius]); %compare to maximum
    longVec1BS = repmat(longVecSMagTmp,probDim,1) .* longVecBSDir; %new long vector
    %at end
    longVecEMagTmp = max([longVecBEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
    longVecEMagTmp = min([longVecEMagTmp;maxSearchRadius]); %compare to maximum
    longVec1BE = repmat(longVecEMagTmp,probDim,1) .* longVecBEDir; %new long vector

    %make sure that "short vectors" at track starts are within
    %allowed range
    shortVecSMagTmp = max([shortVecSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
    shortVecSMagTmp = min([shortVecSMagTmp;maxSearchRadius]); %compare to maximum
    shortVecS1 = repmat(shortVecSMagTmp,probDim,1) .* shortVecSDir; %new short vector

    %make sure that "short vectors" at track ends are within allowed
    %range
    shortVecEMagTmp = max([shortVecEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
    shortVecEMagTmp = min([shortVecEMagTmp;maxSearchRadius]); %compare to maximum
    shortVecE1 = repmat(shortVecEMagTmp,probDim,1) .* shortVecEDir; %new short vector

    %save values for this track
    longVecSAll(:,:,iTrack) = longVec1FS;
    longVecEAll(:,:,iTrack) = longVec1FE;
    shortVecSAll(:,:,iTrack) = shortVecS1;
    shortVecEAll(:,:,iTrack) = shortVecE1;

    longRedVecSAll(:,:,iTrack) = longVec1BS;
    longRedVecEAll(:,:,iTrack) = longVec1BE;

end %(for iTrack = 1 : numTracks)

%% Gap closing

%find all pairs of ends and starts that can potentially be linked
%determine this by looking at time gaps between ends and starts
%and by looking at the distance between pairs
indxEnd2 = [];
indxStart2 = [];

%go over all frames until the one before last
for iFrame = 1 : numFrames - 1
    
    %find tracks that end in this frame
    endsToConsider = tracksPerFrame(iFrame).ends;
    
    for jFrame = iFrame + 1 : min(iFrame+timeWindow,numFrames)
        
        %find tracks that start in this frame
        startsToConsider = tracksPerFrame(jFrame).starts;
        
        %calculate the distance between ends and starts
        dispMat2 = createDistanceMatrix(coordEnd(endsToConsider,:),...
            coordStart(startsToConsider,:));

        tmpFrame = jFrame-iFrame;
        [indxEnd3,indxStart3] = find(dispMat2 <= (maxSpeed * tmpFrame));

        if size(indxEnd3,1) == 1
            indxEnd3 = indxEnd3';
            indxStart3 = indxStart3';
        end

        %add them to the list of possible pairs
        indxEnd2 = [indxEnd2; endsToConsider(indxEnd3)];
        indxStart2 = [indxStart2; startsToConsider(indxStart3)];
                
    end %(for jFrame = iFrame + 1 : iFrame + timeWindow)
    
end %(for iFrame = 1 : numFrames)

%get total number of pairs
numPairs = length(indxEnd2);

%reserve memory for cost matrix vectors
indx1 = zeros(numPairs,1); %row number in cost matrix
indx2 = zeros(numPairs,1); %column number in cost matrix
cost  = zeros(numPairs,1); %cost value

%put time scaling of linear motion in a vector
% timeScalingLin = ones(timeWindow,1);
timeScalingLin = [(1:timeReachConfL).^linScaling(1) ...
    (timeReachConfL)^linScaling(1) * (2:timeWindow-timeReachConfL+1).^linScaling(2)];


%go over all possible pairs of starts and ends
for iPair = 1 : numPairs
    
    %get indices of starts and ends
    iStart = indxStart2(iPair);
    iEnd = indxEnd2(iPair);
    
    %determine the time gap between them
    timeGap = trackStartTime(iStart) - trackEndTime(iEnd);
    
    %calculate the vector connecting the end of track iEnd to the
    %start of track iStart and compute its magnitude
    dispVec = coordStart(iStart,:) - coordEnd(iEnd,:);
    dispVecMag = norm(dispVec);

    %determine whether the connecting vector is parallel or anti-parallel
    %to the tracks' directions of motion
    parallelToS = (dispVec * xyzVelS(iStart,:,1)') > 0;
    parallelToE = (dispVec * xyzVelE(iEnd,:,1)') > 0;
    
    %determine the search area of track iStart
    if ~parallelToS
        longVecS = longRedVecSAll(:,timeGap,iStart);
    else
        longVecS = longVecSAll(:,timeGap,iStart);
    end
    shortVecS = shortVecSAll(:,timeGap,iStart);
        
    %determine the search area of track iEnd
    if ~parallelToE
        longVecE = longRedVecEAll(:,timeGap,iEnd);
    else
        longVecE = longVecEAll(:,timeGap,iEnd);
    end
    shortVecE = shortVecEAll(:,timeGap,iEnd);

    %calculate the magnitudes of the long and short search vectors
    %of both start and end
    longVecMagS = norm(longVecS);
    shortVecMagS = norm(shortVecS);
    longVecMagE = norm(longVecE);
    shortVecMagE = norm(shortVecE);

    %project the connecting vector onto the long and short vectors
    %of track iStart and take absolute value
    projStartLong = abs(dispVec * longVecS) / longVecMagS;
    projStartShort = abs(dispVec * shortVecS) / shortVecMagS;
    
    %project the connecting vector onto the long and short vectors
    %of track iEnd and take absolute value
    projEndLong = abs(dispVec * longVecE) / longVecMagE;
    projEndShort = abs(dispVec * shortVecE) / shortVecMagE;
                  
    %calculate the cosine of the angle between velocity
    %vectors
    cosAngle = longVecE' * longVecS / (longVecMagE * longVecMagS);

    %calculate the square sine of the angle between velocity vectors
    sin2Angle = 1 - cosAngle^2;

    %calculate the square sine of the angle between each
    %motion direction vector and the center-to-center vector
    sin2AngleE = 1 - (dispVec * longVecE / ...
        (longVecMagE * dispVecMag))^2;
    sin2AngleS = 1 - (dispVec * longVecS / ...
        (longVecMagS * dispVecMag))^2;

    %check whether 
    %(1) the end of track iEnd is within the search
    %rectangle of the start of track iStart,
    %(2) the start of track iStart is within the search
    %rectangle of the end of track iEnd,
    %(3) the angle between the two directions of motion
    %is within acceptable bounds, and 
    %(4) the angle between directions of motion and vector
    %connecting end and start is within acceptable bounds
    possibleLink = ...
        ((projEndLong <= longVecMagE) && ...
        (projEndShort <= shortVecMagE)) && ...
        ...
        ((projStartLong <= longVecMagS) && ...
        (projStartShort <= shortVecMagS)) && ...
        ...
        (sin2Angle <= sin2AngleMax) && ...
        ...
        ((sin2AngleE <= sin2AngleMaxVD) && ...
        (sin2AngleS <= sin2AngleMaxVD));

    %only allow links between tracks moving in the same direction
   possibleLink = possibleLink && (cosAngle >= 0);
    
    %if this is a possible link ...
    if possibleLink
        
        %calculate the average displacement for the two tracks combined
        meanDispTrack1 = trackMeanDispMag(iStart);
        meanDispTrack1(isnan(meanDispTrack1)) = meanDispAllTracks;
        meanDispTrack2 = trackMeanDispMag(iEnd);
        meanDispTrack2(isnan(meanDispTrack2)) = meanDispAllTracks;
        meanDisp2Tracks = mean([meanDispTrack1 meanDispTrack2]);
        
        %calculate the cost of linking
        dispVecMag2 = dispVecMag ^ 2;
        cost12 = dispVecMag2 * (1 + mean([sin2Angle sin2AngleE sin2AngleS])) ...
            / (timeScalingLin(timeGap) * meanDisp2Tracks)^2;
        
        %if the lifetime consideration does not make this link impossible
        if isfinite(cost12)
            
            %penalize cost for gap length considerations
            cost12 = cost12 * gapPenalty^(timeGap-1);
            
            %add this cost to the list of costs
            cost(iPair) = cost12;
            
            %specify the location of this pair in the cost matrix
            indx1(iPair) = iEnd; %row number
            indx2(iPair) = iStart; %column number
            
        end
        
    end %(if possibleLink)
    
end %(for iPair = 1 : numPairs)

%keep only pairs that turned out to be possible
possiblePairs = find(indx1 ~= 0);
indx1 = indx1(possiblePairs);
indx2 = indx2(possiblePairs);
cost  = cost(possiblePairs);


%% Append cost matrix to allow births and deaths ...

%determine the cost of birth and death
tmp = (costMat~=0);
numPotAssignRow = full(sum(tmp,2));
numPotAssignCol = full(sum(tmp)');
numPotAssignColAll = sum(numPotAssignCol);
numPotAssignRowAll = sum(numPotAssignRow);
numPartCol = length(numPotAssignCol) * 2;
extraCol = (numPotAssignColAll-numPartCol)/numPotAssignColAll;
numPartRow = length(numPotAssignRow) * 2;
extraRow = (numPotAssignRowAll-numPartRow)/numPotAssignRowAll;
prctile2use = min(100, 100 - mean([extraRow extraCol])*100);

costBD = 1.05*prctile(cost(:),prctile2use);

%get the cost for the lower right block
% costLR = min(min(min(costMat))-1,-1);
costLR = costBD;

%create cost matrix that allows for births and deaths
% costMat = [costMat ... %costs for links (gap closing + merge/split)
%     spdiags([costBD*ones(numTracks,1); altCostSplit],0,numEndSplit,numEndSplit); ... %costs for death
%     spdiags([costBD*ones(numTracks,1); altCostMerge],0,numStartMerge,numStartMerge) ...  %costs for birth
%     sparse(indx2,indx1,costLR*ones(length(indx1),1),numStartMerge,numEndSplit)]; %dummy costs to complete the cost matrix

costMat = [costMat ... %costs for links (gap closing + merge/split)
    spdiags(costBD*ones(numTracks+numSplit,1),0,numEndSplit,numEndSplit); ... %costs for death
    spdiags(costBD*ones(numTracks+numMerge,1),0,numStartMerge,numStartMerge) ...  %costs for birth
    sparse(indx2,indx1,costLR*ones(length(indx1),1),numStartMerge,numEndSplit)]; %dummy costs to complete the cost matrix

%determine the nonlinkMarker
nonlinkMarker = min(floor(full(min(min(costMat))))-5,-5);
    
    
end

