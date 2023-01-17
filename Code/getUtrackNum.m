function numUtracks = getUtrackNum(utrackPackPaths)
%GETUTRACKNUM: look for utrack results in folder and count number of compound tracks
%
%SYNOPSIS: numUtracks = getUtrackNum(utrackPackPaths)
%
%INPUT: utrackPackPaths : (nx1 cell array) directories to utrack results
%
%OUTPUT:     numUtracks : (nx1 cell array) number of compound tracks in each
%                         result
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

%% Input control & variable initiation
numMov = size(utrackPackPaths,1);
fs = filesep;

%% Output
numUtracksTmp = cell(numMov,1);

%% Loop through each movie utrack result and count number of compound track
for i = 1:numMov
    % STEP 0: if one of the paths are empty, move on to the next listed path
    if isempty(utrackPackPaths{i,1})
        continue
    end
    
    % STEP 1: Search for and Load utrack result
    
    % load uTrack output in "TrackingPackage/tracks" folder
    try
        load([utrackPackPaths{i,1} fs 'tracks' fs 'Channel_1_tracking_result.mat'], 'tracksFinal');
    catch
        error('tracksFinal DOES NOT LOAD!!!!!!!!')
    end
    
    % STEP 2: Count
    numUtracksTmp{i,1} = size(tracksFinal,1);
    
end

numUtracks = cell2mat(numUtracksTmp);
end