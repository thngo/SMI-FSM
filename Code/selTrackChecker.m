function [passTrackLogic, numPassTrack, statPassTrack, lftPassTrack] = selTrackChecker(decompTrack, checkFr, minLft)
%PASSTRACKCHECKER checks for decompound tracks that pass through pre-defined frame and satisfy minimum lifetime
%
%SYNOPSIS: [passTrack, passTrackId] = selTrackChecker(decompTrack, frameCheck, minLft)
%
%INPUT:     decompTrack: (n x 1 struct) ouput of decompoundCompTrack.m
%
%              checkFr : (1 x m vector of integer) the different frames to
%                        check whether decompTrack pass through them
%
%               minLft : (integer) a track is required to have > minLft.
%
%OUTPUT: passTrackLogic: (n x 1 logic vector) 0 = indicate the ith track
%                        did not pass the conditions that was input; 1 =
%                        indicate that the ith track passed the conditions
%                        that were input.
%
%          numPassTrack: (integer) number of tracks that passed the
%                        conditions that were input.
%
%         statPassTrack: (3 x 1 vector) statistics [min max median] length 
%                        of the tracks that satisfy the input conditions. 
%               Note that "tracks that pass through frame X" means that
%               these tracks started before and after frame X. If a track
%               starts at frame X or ends at frame X, it would not qualify.
%
%          lftPassTrack: (n x 1 vector) array of lifetimes of tracks that 
%                        passed the input conditions.
%
%Tra Ngo, July 2022
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

%% Input &
if nargin < 1
    error('Expect at least 1 input.')
end

if ~exist('checkFr','var')
    checkFr = [];
end

if ~exist('minLft','var')
    minLft = [];
end

numCheck = length(checkFr);
numTrack = size(decompTrack,1);

%% Initialize output:
passTrackLogic = [];
statPassTrack = nan(3,1);

%% Get track SEL
trackSEL = getTrackSEL(decompTrack);

% CONDITION 1: Starts before border frame & ends after border frame
if ~isempty(checkFr)
    checkFrLogic = false(numTrack,numCheck);
    for i = 1:numCheck
        checkFrLogic(:,i) = trackSEL(:,1) < checkFr(i) & trackSEL(:,2) > checkFr(i);
    end
    checkAllFrLogic = any(checkFrLogic,2); % satisfy condition 1 for checkFr(i) or checkFr(2) or so on.
else
    checkAllFrLogic = true(numTrack,1); % take all frame if there was no checking
end

% CONDITION 2: Satisfy minimum lifetime:
if ~isempty(minLft)
    lftLogic = trackSEL(:,3) > minLft;
else
    lftLogic = true(numTrack,1); % take all frame if there was no checking
end


% Logical array of tracks that satisfy all conditions:
passTrackLogic = checkAllFrLogic & lftLogic;

numPassTrack = sum(passTrackLogic);

% Get lifetime vector:
lftPassTrack = trackSEL(passTrackLogic,3);

% Get statistics of the tracks that satisfy the conditions.
statPassTrack(1,1) = min(lftPassTrack);
statPassTrack(2,1) = max(lftPassTrack);
statPassTrack(3,1) = median(lftPassTrack);

end