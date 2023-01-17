function tracksCrop = cropCompTracksFR(tracks0,frameRange)
%cropCompTracksFR crops tracks to a particular frame range
%
%SYNOPSIS tracksCrop = cropCompTracksFR(tracks0,frameRange)
%
%INPUT  tracks0      : Original tracks, in form of output of
%                      trackCloseGapsKalman.
%       frameRange   : Row vector with two entries indicating time range to
%                      retain.
%
%OUTPUT tracksCrop   : Same as input tracks, just cropped in time.
%
%REMARKS Code still needs to handle case when a compound track gets broken
%        into multiple non-interacting tracks.
%
%REMARKS2: Aggregation state is not handled w.r.t. sparse form. -TN20200920
%
%Khuloud Jaqaman, December 2014
%
% Copyright (C) 2022, Jaqaman Lab - UTSouthwestern 
%
% This file is part of SMI-FSM.
% 
% SMI-FSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SMI-FSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SMI-FSM.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%% Input

if nargin < 2
    error('cropCompTracksFR: Incorrect number of input arguments!')
end


%check if tracks are in sparse form
sparseForm = issparse(tracks0(1).tracksFeatIndxCG);

%check if input contain aggregation state (which will need to be processed
%if exist)
aggregField = isfield(tracks0,'aggregState');

%% Cropping

%get number of tracks
numTracks = length(tracks0);

%go over all tracks and crop frames
tracksCrop = tracks0;
for iTrack = 1 : numTracks
    
    
    %convert the current track's information to matrix format
    [trackedFeatureInfo,trackedFeatureIndx, ~, ~,aggregStateMat] = convStruct2MatIgnoreMS(tracks0(iTrack));
    
    %determine which frames to keep by comparing the last time point of this
    %compound track from the number of column vs. the input frameRange(2) -- TN20190924
    tpKeep = frameRange(1) : min(size(trackedFeatureIndx,2), frameRange(2));
    
    %thus determine which columns to keep in tracksCoordAmpCG
    colKeep = (repmat(tpKeep,8,1)-1)*8 + repmat((1:8)',1,length(tpKeep));
    colKeep = colKeep(:);
    
    %keep only the time points of interest
    trackedFeatureIndx = trackedFeatureIndx(:,tpKeep);
    trackedFeatureInfo = trackedFeatureInfo(:,colKeep);
    if aggregField
        aggregStateMat = aggregStateMat(:,tpKeep);
    end
    
    %convert zeros to NaNs if original tracks were sparse
    if sparseForm
        trackedFeatureInfo(trackedFeatureInfo==0) = NaN;
    end
    
    %get each track segment's new start, end and life time
    segSEL = getTrackSEL(trackedFeatureInfo); % start time, end time, lifetime
    numSeg = size(segSEL,1);
    
    % Remove the padding 0s preceeding the start of an actual track -- TN20190924
    % Remove the run 0s after the ultimate end of an actual compTrack -- TN20190927
    ctStartTime = min(segSEL(:,1)); % start time of compTrack
    ctEndTime = max(segSEL(:,2)); % end time of compTrack
    if any(segSEL(:,1))
        trackedFeatureIndx = trackedFeatureIndx(:,ctStartTime:ctEndTime);
        trackedFeatureInfo = trackedFeatureInfo(:,(8*(ctStartTime-1)+1):(8*ctEndTime));
        if aggregField
            aggregStateMat = aggregStateMat(:,ctStartTime:ctEndTime);
        end
    end
    
    % add the shift back to segSEL so that seqOfEvents follows the original
    % compound track's time frame -- TN20190924
    segSEL = segSEL + frameRange(1) - 1;
    
    %find segments that survive the cropping and those that do not
    indxStay = find(~isnan(segSEL(:,3)));
    indxGone = setdiff((1:numSeg)',indxStay);
    if isempty(indxGone)
        indxGone  = [];
    end
    
    %get the track's original sequence of events
    seqOfEvents = tracks0(iTrack).seqOfEvents;
    
    %"delete" merges and splits happening outside of the frame range
    tmp = ~isnan(seqOfEvents(:,4)) & ( seqOfEvents(:,1)<frameRange(1) | seqOfEvents(:,1)>frameRange(2) );
    seqOfEvents(tmp,4) = NaN;
    
    %assign new start and end times - merges will need an addition of one,
    %done two steps down
    for iSeg = 1 : numSeg
        rowsSeg = find(seqOfEvents(:,3)==iSeg);
        
        %%%%% Modified by Luciana de Oliveira September 2016
        % There are some events that does not have the initial time, so for
        % these cases the vector has only one value and Matlab point an
        % error case, for solve that I changed and put an if condition.
        
        if length(rowsSeg)>1
            seqOfEvents(rowsSeg(1),1) = segSEL(iSeg,1);
            seqOfEvents(rowsSeg(2),1) = segSEL(iSeg,2);
        else
            seqOfEvents(rowsSeg(1),1) = segSEL(iSeg,1);
        end
    end
    
    %replace segments that did not survive the cropping with NaN in both
    %the 3rd and 4th column
    for iSeg = indxGone'
        seqOfEvents(seqOfEvents(:,3)==iSeg,3) = NaN;
        seqOfEvents(seqOfEvents(:,4)==iSeg,4) = NaN;
    end
    
    %for surviving merges, add one to end time to follow convention
    rowsSeg = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)));
    seqOfEvents(rowsSeg,1) = seqOfEvents(rowsSeg,1) + 1;
    
    %keep only rows that belong to surviving segments
    seqOfEvents = seqOfEvents(~isnan(seqOfEvents(:,3)),:);
    trackedFeatureInfo = trackedFeatureInfo(indxStay,:);
    trackedFeatureIndx = trackedFeatureIndx(indxStay,:);
    if aggregField
        aggregStateMat = aggregStateMat(indxStay,:);
    end
    
    %renumber segments to reflect new trackedFeatureInfo and
    %trackedFeatureIndx
    for iStay = 1 : length(indxStay)
        iSeg = indxStay(iStay);
        seqOfEvents(seqOfEvents(:,3)==iSeg,3) = iStay;
        seqOfEvents(seqOfEvents(:,4)==iSeg,4) = iStay;
    end
    
    % rearrage the sequence of events by starting frame - TN20190930
    % Remark: mind that SM (or track) appearance must preceed disappearance
    % and currently "sort" function seems to not disturb this ordering.
    % However, we have no built-ins to make sure that this happens.
    [~,idx] = sort(seqOfEvents(:,1));
    seqOfEvents_final = seqOfEvents(idx,:);
    
    %convert to sparse if input was sparse
    if sparseForm
        trackedFeatureIndx = sparse(trackedFeatureIndx);
        trackedFeatureInfo(isnan(trackedFeatureInfo)) = 0;
        trackedFeatureInfo = sparse(trackedFeatureInfo);
    end
    
    %store the cropped compound track
    tracksCrop(iTrack).tracksFeatIndxCG = trackedFeatureIndx;
    tracksCrop(iTrack).tracksCoordAmpCG = trackedFeatureInfo;
    tracksCrop(iTrack).seqOfEvents = seqOfEvents_final;
    if aggregField
        tracksCrop(iTrack).aggregState = aggregStateMat;
    end
    
end

%% Modification LRO 2019/04/11- to reagroup compTracks after resampling
%reagroup compTracks
tracksCrop = groupSegmentsComptracksAfterReformation(tracksCrop);

%% remove the empty tracks
flagKeep = true(length(tracksCrop),1);
for iTrack = 1:length(tracksCrop)
    flagKeep(iTrack) = ~isempty(tracksCrop(iTrack).tracksFeatIndxCG);
end
tracksCrop = tracksCrop(flagKeep);

%% ~~~ the end ~~~