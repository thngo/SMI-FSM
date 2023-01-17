function mTrackChopped = simpleTrackChop(simpleTrack, m)
% SIMPLETRACKCHOP : chops the simpleTrack into m tracklets and return the chopped tracklets.
%
% SYNOPSIS: mTrackChopped = simpleTrackChop(simpleTrack, m)
%
% INPUT:    simpleTrack  : (structure) contains 1 row (1 track) with fields
%                           .tracksFeatIndxCG
%                           .tracksCoordAmpCG
%                           .seqOfEvents
%                           .aggregState (optional, will be handled if
%                           input includes this field)
%                   m    : (integer > 0) number of resulting tracklets.
%                           Default: m = 1
%
% OUTPUT:  mTrackChopped : (structure) each row contains the tracklets
%                           that simpleTrack was chopped into, with fields
%                           .tracksFeatIndxCG: Leading & trailing 0s trimmed.
%                           .tracksCoordAmpCG: Leading & trailing NaNs trimmed.
%                           .seqOfEvents: all events are true event ("NaN");
%                                         all tracks are numbered "1";
%                                         start & end frames comes from simpleTrack.
%                           .aggregState (optional, will be handled if
%                           input includes this field)
%
% Tra Ngo, Mar 2020
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

%% Input control
if nargin < 2 || m < 1
    m = 1;
    disp('Input track is by default not being chopped. Please input valid number of chop intervals!!!')
end

if isempty(simpleTrack)
    error('Input track is empty!!!')
end

%check if input contain aggregation state (which will need to be processed
%if exist)
aggregField = isfield(simpleTrack,'aggregState');

%% Output initialization
% For small m (<100), using for-loop to initialize structure mTrackChopped is faster than using repmat.
% mTrackChopped2 = repmat(struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],'seqOfEvents',[]),m,1);

% for i = 1:m
%     mTrackChopped(i,1).tracksFeatIndxCG = [];
%     mTrackChopped(i,1).tracksCoordAmpCG = [];
%     mTrackChopped(i,1).seqOfEvents = [];
%     if aggregField
%         mTrackChopped(i,1).aggregState = [];
%     end
% end

mTrackChopped(1,1).tracksFeatIndxCG = [];
mTrackChopped(1,1).tracksCoordAmpCG = [];
mTrackChopped(1,1).seqOfEvents = [];
if aggregField
    mTrackChopped(1,1).aggregState = [];
end

%% Determine chopped tracklets' size
divTrackSize = size(simpleTrack.tracksFeatIndxCG,2)/m;
trackLength = size(simpleTrack.tracksFeatIndxCG,2);

rowCorrector = 0;
for iAdd = 1:m
    
    % Determine start frame and end frame (trimming leading and trailing 0s
    % and NaNs of each tracklet)
    startFrame = max(1,ceil((iAdd-1)*divTrackSize));
    endFrame = min(ceil(divTrackSize*iAdd),trackLength);
    if sum(simpleTrack.tracksFeatIndxCG(1,startFrame:endFrame)) == 0
        % if tracks is only 0s, skip
        rowCorrector = rowCorrector +1;
    else
        while simpleTrack.tracksFeatIndxCG(1,(startFrame)) == 0
            startFrame = startFrame + 1;
        end
        
        while simpleTrack.tracksFeatIndxCG(1,(endFrame)) == 0
            endFrame = endFrame - 1;
        end
        
        
        % tracksFeatIndxCG & tracksCoordAmpCG:
        mTrackChopped(iAdd-rowCorrector).tracksFeatIndxCG = ...
            simpleTrack.tracksFeatIndxCG(1, (startFrame):(endFrame));
        mTrackChopped(iAdd-rowCorrector).tracksCoordAmpCG = ...
            simpleTrack.tracksCoordAmpCG(1,((startFrame-1)*8+1):(endFrame*8));
        
        % TN 20200925: add handling of aggregState
        if aggregField
            mTrackChopped(iAdd-rowCorrector).aggregState = ...
                simpleTrack.aggregState(1, (startFrame):(endFrame));
        end
        
        
        % create new .seqOfEvents (soe)
        soeStartFrame = simpleTrack.seqOfEvents(1,1) + startFrame - 1;
        soeEndFrame = simpleTrack.seqOfEvents(1,1) + endFrame - 1;
        tmpSeqOfEvents = [soeStartFrame, 1, 1, NaN; soeEndFrame, 2, 1, NaN];
        
        mTrackChopped(iAdd-rowCorrector).seqOfEvents = tmpSeqOfEvents;
    end
end

end