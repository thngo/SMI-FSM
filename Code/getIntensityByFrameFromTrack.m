function [coltInt] = getIntensityByFrameFromTrack(tracksClean, endFrame, beginFrame)
%GETINTENSITYBYFRAMEFROMTRACK User inputs X number of frame and function outputs column vector of intensities of detections from first X frames of all tracks
%
%SYNOPSIS:  [coltInt] = getIntensityByFrameFromTrack(tracksClean, endFrame, beginFrame)
%
%INPUT:         tracksClean: (structure) compound track after clean up with
%                            removeSimultaneousMergeSplit,
%                            reformatSplitsMergesKeepLongerLivedSegment,
%                            removeSplitMergeArtifactsChronological(,[])
%                            each track can include multiple interacting
%                            segments; includes following fields:
%                               .tracksFeatIndxCG
%                               .tracksCoordAmpCG
%                               .seqOfEvents
%
%                   endFrame  : (optional, integer) frame X at which 
%                            collection of intensities starts (default frame No.1)
%
%                   beginFrame: (optional, integer) frame Y at which
%                            collection of intensities ends (default frame No. 5)
%
%OUTPUT:        coltInt: column vector of intensities collected from first
%                        X frames included in tracksClean.
%
%Tra Ngo, March 2020
%% Initialize output:
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
coltInt = [];

%% Check input
if isempty(tracksClean) || ~isfield(tracksClean,"tracksFeatIndxCG") || ~isfield(tracksClean,"tracksCoordAmpCG") || ~isfield(tracksClean,"seqOfEvents")
    error("Input compound track format invalid!");
end

if nargin < 2
    beginFrame = 1;
    endFrame = 5; % take the first 5 frames for intensity fit calculation
elseif nargin < 3
    beginFrame = 1;
elseif nargin == 3 && endFrame < beginFrame
    disp("!!!!!!!!! Input end frame (no.2) is before begin frame (no.3) !!!!!!!!!");
end

%% Collect intensities:
% Method 1:
% Compound tracks are decompounded so that each row contains only 1 tracklet
decompTrack = decompoundCompTracks(tracksClean);
trackSEL = getTrackSEL(decompTrack);

% tracks starts before set numFr and ends at or after set numFr
indxTrack1st5 = find(trackSEL(:,1) <= endFrame & trackSEL(:,2) >= beginFrame);

for iTrk = indxTrack1st5'
    % count starting at beginFrame or first frame of track if first frame is
    % after beginFrame
    beginFrameCurr = max(1, beginFrame-trackSEL(iTrk,1)+1);
    % count until endFrame or last frame of track if last frame is before
    % endFrame
    endingFrameCurr = min(trackSEL(iTrk,3), endFrame-trackSEL(iTrk,1)+1); 
    
    for iFr = beginFrameCurr:endingFrameCurr
        tmpColtInt = decompTrack(iTrk).tracksCoordAmpCG(:,(iFr-1)*8+4); % intensity information is in the 4th other col.
        tmpColtInt = tmpColtInt(~isnan(tmpColtInt),:);
        coltInt = vertcat(coltInt, tmpColtInt);
    end
end


% % % % Method 2: seems to take double the time of method 1
% % % tmpTrackMat = convStruct2MatIgnoreMS(tracksClean);
% % % coltInt2 = [];
% % % for iFr = 1:numFr
% % %     tmpColtInt = tmpTrackMat(:,(iFr-1)*8+4);
% % %     tmpColtInt = tmpColtInt(~isnan(tmpColtInt),:);
% % %     coltInt2 = vertcat(coltInt2, tmpColtInt);
% % % end

end

