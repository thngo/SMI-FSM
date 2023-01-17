function [trackCntOut] = getFsmTrackNum(qfsmPackPaths, frameRange)
%GETFSMTRACKNUM: output number of speckle track per movie from MPM of QFSM given directory of movie where QFSM was performed
%
%SYNOPSIS: [trackCntOut] = getFsmTrackNum(qfsmPackPaths, frameRange)
%
%INPUT : qfsmPackPaths : (nx1 cell) locations to n actin movies
%                        Note that paths in input must be appropriate to OS
%                        calling this function (i.e. Linux paths on Linux
%                        OS, Windows paths on Windows OS)
%                        i.e. {1,1} = {'/project/biophysics/jaqaman_lab/vegf_tsp1/adasgu/FSM_SMI_ABD_ABDRA/20200113_PD80_LP70_Transfex_point25ug_CD36_C2_D2_day2_halo_3nMJF549_TIMEmngActlowP14/CD36/Actin/m-01-02'}
%
%           frameRange : (2x1 vector, optional) start frame & end frame where 
%                        number of speckle tracks are counted; default =
%                        the whole range of available frames from start:end
%
%OUTPUT: trackCntOut   : (nx1 cell) number of speckle track within each
%                         input movie.
%
% Tra Ngo, July 2022
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

%% Input and initialization:
numMovs = length(qfsmPackPaths);
fs = filesep;

if ~exist('frameRange','var') || isempty(frameRange)
    frameRange = [];
elseif size(frameRange,1) ~= 1 && size(frameRange,2) ~= 2
    error('Expected frameRange input to be 1x2 vector of integer.')
end
    
for i = 1 : numMovs
    
trackCnt = 0;
iEntry = 1;
   
    %skip an element in the cell if path does not exist.
    if isempty(qfsmPackPaths{i,1})
        continue
    end
    
    %% Load MPM
    % Load speckle tracks (MPM)
    load([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'speckleTracks' fs '_tracks.mat'], 'MPM');
    
    if ~isempty(frameRange) % if user input a range of frame for trimming MPM
        startFr = frameRange(1,1)*2 - 1;
        endFr = frameRange(1,2)*2;
        MPM = MPM(:,(startFr:endFr));
    end
    
    [~,lifetimeMPM] = retainMPMminLifetime(MPM, 2);
    % Gather only the x-pos on lifetimeMPM
    lifetimeMPM_x = lifetimeMPM(:,1:2:end);
    %notGhostMPM_x_logic = lifetimeMPM_x > 1;
    %tmp = diff(notGhostMPM_x_logic);
    %sum(lifetimeMPM_x,2);
    
    %% count speckle tracks
    
    % Loop through each entry in MPM, check the track's lifetime, then skip to next
    % track
    lifetimeMPM_x_0 = [lifetimeMPM_x 0*ones(size(lifetimeMPM_x,1),1)];
    lifetimeMPM_x_0_T = lifetimeMPM_x_0'; % transpose so that we can traverse through the rows of MPM_x
    
    while iEntry <= length(lifetimeMPM_x_0_T(:))
        if lifetimeMPM_x_0_T(iEntry) > 1
            trackCnt = trackCnt + 1;
            iTmp{i}{iEntry} = lifetimeMPM_x_0_T(iEntry);
            iEntry = iEntry + lifetimeMPM_x_0_T(iEntry)+1;
        end
        iEntry = iEntry + 1;
    end % while
    
    trackCntOut{i,1} = trackCnt;
end % for i = 1:numMovs
end