function [shuffledTassm] = shuffleTassmSpkl(tassm,flagShuffleBtwGr)
%SHUFFLETASSMSPKL  : takes tassm and shuffle the rows of speckle properties observations within trajectory groups
%
%SYNOPSIS: [shuffledTassm] = shuffleTassmSpkl(tassm)
%
%INPUT:            tassm : cell array containing :
%                       tassm{:,1}{1,1}  : speckle properties
%                       tassm{:,1}{1,2}  : SM properties
%               with rows be different observations and col be different
%               properties. Observation of speckle is matched with
%               observation of SM in similar row.
%       flagShuffleBtwGr : (logic)
%                 true  = shuffle observations within and between groups.
%                 false = shuffle observations only within groups (default)
%
%OUTPUT: similar to tassm but observation of speckle is randomized (NOT
%        matched with) observation of SM in similar row.
%
% Huong Tra Ngo, Nov 2019
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

if ~exist("flagShuffleBtwGr ",'var') || isempty(flagShuffleBtwGr)
    flagShuffleBtwGr = false;
end

shuffledTassm = tassm;

%% Looping through each cell in tassm array
for iMov = 1:size(tassm,1)
    
    % skip cell in array if cell is empty.
    if isempty(tassm{iMov,1})
        continue
    end
    
    % Shuffle within trajClass group
    for iDiff = 0 : 3
        iTrajInd = find(tassm{iMov,1}{1,2}(:,end) == iDiff); % get original index of observation in this trajClass
        nTrajObs = size(iTrajInd,1); % get number of observation in this trajClass
        randTrajVect = randperm(nTrajObs,nTrajObs); % create randomized index vector
        randTrajInd = iTrajInd(randTrajVect,1); % shuffle original index of observation in this trajClass
        
        % replace speckle properties in proper matched order with shuffled (randomized) order
        shuffledTassm{iMov,1}{1,1}(iTrajInd,:) = tassm{iMov,1}{1,1}(randTrajInd,:);
    end
    
    % Shuffle between different trajClass group
    if flagShuffleBtwGr
        nSpklObservation = size(tassm{iMov,1}{1,1},1); % get number of observation
        randVect = randperm(nSpklObservation,nSpklObservation); % create randomized index vector
        
        % replace speckle properties in proper matched order with shuffled (randomized) order
        shuffledTassm{iMov,1}{1,1} = tassm{iMov,1}{1,1}(randVect,:);
    end
end % (iMov = 1:size(tassm,1))

end