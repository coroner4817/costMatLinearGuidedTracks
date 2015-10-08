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

% get user-set parameters
tMax = gapCloseParam.timeWindow;

%find the number of tracks to be linked and the number of frames in the movie
[nTracks,nFrames] = size(trackedFeatInfo);
nFrames = nFrames / 8;

%list the tracks that start and end in each frame
tracksPerFrame = repmat(struct('starts',[],'ends',[]),nFrames,1);
for iFrame = 1 : nFrames
    tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
    tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
end

%% Gap closing

% extract feature positions and velocity components
px=trackedFeatInfo(:,1:8:end);
py=trackedFeatInfo(:,2:8:end);
vx=diff(px,1,2);
vy=diff(py,1,2);

% TRACK STARTS
trackStartPxyVxy = zeros(nTracks,4);
% x and y coordinates of the track's first point
trackStartPxyVxy(:,1)=cell2mat(arrayfun(@(i) px(i,find(~isnan(px(i,:)),1,'first')),[1:nTracks]','UniformOutput',0));
trackStartPxyVxy(:,2)=cell2mat(arrayfun(@(i) py(i,find(~isnan(py(i,:)),1,'first')),[1:nTracks]','UniformOutput',0));
% average of first three velocity vectors (made from last 4 points
% on track, if that many exist), x and y components
trackStartPxyVxy(:,3)=cell2mat(arrayfun(@(i) mean(vx(i,find(~isnan(vx(i,:)),3,'first'))),[1:nTracks]','UniformOutput',0));
trackStartPxyVxy(:,4)=cell2mat(arrayfun(@(i) mean(vy(i,find(~isnan(vy(i,:)),3,'first'))),[1:nTracks]','UniformOutput',0));

% TRACK ENDS
trackEndPxyVxy = zeros(nTracks,4);
% x and y coordinates of the track's last point
trackEndPxyVxy(:,1)=cell2mat(arrayfun(@(i) px(i,find(~isnan(px(i,:)),1,'last')),[1:nTracks]','UniformOutput',0));
trackEndPxyVxy(:,2)=cell2mat(arrayfun(@(i) py(i,find(~isnan(py(i,:)),1,'last')),[1:nTracks]','UniformOutput',0));
% average of last three velocity vectors (made from last 4 points
% on track, if that many exist), x and y components
trackEndPxyVxy(:,3)=cell2mat(arrayfun(@(i) mean(vx(i,find(~isnan(vx(i,:)),3,'last'))),[1:nTracks]','UniformOutput',0));
trackEndPxyVxy(:,4)=cell2mat(arrayfun(@(i) mean(vy(i,find(~isnan(vy(i,:)),3,'last'))),[1:nTracks]','UniformOutput',0));

% get velocity components for each track from kalman filter (very similar to trackEndVxy)
xyzVel=cell2mat(arrayfun(@(iTrack) kalmanFilterInfo(trackEndTime(iTrack))...
    .stateVec(trackedFeatIndx(iTrack,trackEndTime(iTrack)),probDim+1:2*probDim),...
    [1:nTracks]','UniformOutput',0));

trackEndSpeed=sqrt(sum(xyzVel.^2,2));
vMax=prctile(trackEndSpeed,95);
vMed=median(trackEndSpeed);

% get start and end frames for each track
sFrameAll=zeros(nTracks,1);
eFrameAll=zeros(nTracks,1);
for iFrame=1:nFrames
    sFrameAll(tracksPerFrame(iFrame).starts)=iFrame;
    eFrameAll(tracksPerFrame(iFrame).ends)=iFrame;
end

% initialize matrices for pair indices and cost components
indx1 = zeros(10*nTracks,1);
indx2 = zeros(10*nTracks,1);
costComponents  = zeros(10*nTracks,5);

linkCount = 1;
    
    
end

