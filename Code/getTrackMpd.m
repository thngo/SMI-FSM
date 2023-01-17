function [diffRadOut, detectPairOfMpd, mpdOriginal] = getTrackMpd(tracksIn)
%GETTRACKMPD: calculate maximum pairwise distance of each detections in single-molecule tracks
%
%SYNOPSIS: [diffRadOut, detectPairOfMpd, mpdOriginal] = getTrackMpd(tracksIn)
%
%INPUT:     tracksIn   : (struct) input tracks, has at least following fields:
%                           .tracksFeatIndxCG
%                           .tracksCoordAmpCG
%                           .seqOfEvents
%                        tracksIn is expected to be decompounded. See
%                        details of decomponded track at: decompoundCompTracks.m
%
%OUTPUT:    diffRadOut : (vector) each row correspond to diffusion radius
%                        of each track in tracksIn (= NaN if track consists
%                        of only 1 detection), estimated by:
%               STEP 1: Get the MPDs of each tracklets
%               STEP 2: Get the max of the MPDs.
%
%       detectPairOfMpd: (n x 2 cell array) coordinates of 2 localization
%                        that resulted in reported MPD
%
%           mpdOriginal: (n x 1 cell array) pairwise distance of every pair of
%                        localization in a track (the maximum is reported in
%                        diffRadOut)
%
%Tra Ngo, Apr 2021
% Dev note: this function used to be a sub-function within file basicTransientTrackChopBySpan.m 
% to get track maximum pair distance of each detections.
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

%% Input & initiate variables
%get number of track segments to be analyzed
numTrackSegments = length(tracksIn);

%get track segment start times, end times and life times
tracks = convStruct2MatIgnoreMS(tracksIn);
trackSEL = getTrackSEL(tracks);

%% Output
diffRadOut = nan(numTrackSegments,1);
detectPairOfMpd = cell(numTrackSegments,2);
mpdOriginal = cell(numTrackSegments,1);

%% Go over all analyzable track segments
for iTrack = 1:numTrackSegments
    
    % Get the current track's S E L (start, end, length)
    trackCurrLength = trackSEL(iTrack,3);
    trackCurrStart = trackSEL(iTrack,1);
    trackCurrEnd = trackSEL(iTrack,2);
    
    if trackCurrLength > 1
        
        %initialize detection vector
        detect = NaN(1,2);
        
        coordAmpCG = tracks(iTrack,8*(trackCurrStart-1)+1:8*trackCurrEnd);
        xCoord = coordAmpCG(1:8:end);
        yCoord = coordAmpCG(2:8:end);
        X = [xCoord',yCoord'];
        
        % STEP 1: Get the MPDs of each tracklets if track's length is larger than 1
        D = pdist(X,'euclidean');
        
        % STEP 2: take the max MPD
        diffRadOut(iTrack) = max(D);
        
        % get the pair of detection that is responsible for MPD
        sqrTmp = squareform(D);
        [detect(1), detect(2)] = find(sqrTmp  == diffRadOut(iTrack),1,'first');
        detectCoord = [xCoord(detect(1)) yCoord(detect(1)); xCoord(detect(2)) yCoord(detect(2)) ];
        
        detectPairOfMpd{iTrack,1} = detect; % detection number A & number B
        detectPairOfMpd{iTrack,2} = detectCoord; % coordinate of detection number A & number B
        
        % Save maxDisplacement for output
        mpdOriginal{iTrack} = D;
        
    else
        %if track consists of only 1 detection, output the MPD & related info as NaN.
        diffRadOut(iTrack) = NaN;
        detectPairOfMpd{iTrack,1} = nan(1,2);
        detectPairOfMpd{iTrack,2} = nan(2,2);
        mpdOriginal{iTrack} = NaN;
        
    end % (if trackPartLength > 1)
end
end