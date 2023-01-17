function utrackPackPaths = getUtrackPackPaths(qfsmPackPaths)
%GETUTRACKPACKPATHS returns the corresponding paths to results of u-track analysis of single-molecules movies for each FSM movie path
%
%SYNOPSIS: utrackPackPaths = getUtrackPackPaths(qfsmPackPaths)
%
%INPUT      qfsmPackPaths: (n x 1 cell array) containing paths to the movie 
%                          folder containing QFSM results
%
%OUTPUT   utrackPackPaths: (n x 1 cell array) containing paths to the "TrackingPackage"
%                          folder containing u-track analysis results for each SMI 
%                          movie corresponding to the FSM movie in the qfsmPackPaths input.
%
%ASSUMPTION: this function assumes a very folder specific structure between
%results of u-track and of QFSM softwares. Same folder must contain the
%receptorFolder and actinFolder, within each must contain sub-folder with
%identical names (e.g. m-01-01, m-01-02, etc.), within each must contain
%"TrackingPackage" and "QFSMPackage" sub-folders respectively.
%
% Tra H. Ngo, February 2019
% Tra H. Ngo, Jun 2019, modified to consider multiple ways of writing "actin" and "m-0x-0x" folder.
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

%% List all possible folder names:
receptorFolder = ["CD36","Cd36","cd36",...
    "ABD_C2","C2_ABD","ABDRA_D2","ABD",...
    "receptors", "receptor", "Receptor", "Receptors", ...
    "singleMolecule", "SingleMolecule", "SMI", "smi",...
    "singleMolecules", "SingleMolecules"]; % name of receptors/SM/surface proteins folders
actinFolder = ["actin","Actin","speckle","Speckle","speckles","Speckles"];

%% Initialize variables and outputs:
utrackPackPaths_cat{1,1} = [];
fs = filesep;

%% Locate folders containing u-track result 
% The following steps assume a very specific folder structure between 
% results of u-track and of QFSM softwares.

for iPwd = 1:length(qfsmPackPaths)
    % navigate to folder that contains the two folders where FSM and u-track were performed    
    cd(qfsmPackPaths{iPwd,1})
    cd ..
    cd ..
    
    % get movie name "-0x-0x" of movie in each row of qfsmPackPaths
    movName = extractBetween(qfsmPackPaths{iPwd,1},strcat(fs,actinFolder(isfolder(actinFolder)),fs,'m'),length(qfsmPackPaths{iPwd,1}));

        %strcat('/',actinFolder(isfolder(actinFolder)) replaces '/Actin/m'
        %in case the 'Actin' folder is not capitalized and written as
        %'actin' instead. Old code as below:
        %movName = extractBetween(qfsmPackPaths{iPwd,1},'/Actin/m',length(qfsmPackPaths{iPwd,1}));
    
    % Navigate inside single-molecule folder
    cd( char(receptorFolder(isfolder(receptorFolder))) )
    baseDir =  (cd);

    % get all names inside single-molecule folder that starts with "m-0x-0x"
    mergeDir = dir(fullfile(baseDir,strcat('m',movName{1,1}(1:6),'*')));
    % get path of folder with name "m-0x-0x"
    mergeDir = fullfile(baseDir,{mergeDir([mergeDir.isdir]).name});
       
    utrackPackPaths_cat{iPwd,1} = string(strcat(mergeDir,filesep,'TrackingPackage'));
    switch length(utrackPackPaths_cat{iPwd,1}) 
        
        case 0 % if there is no "m-0x-0x" folder, we cannot go to "TrackingPackage" folder
            error('Error. \n There does not seem to exist a movie folder corresponding to path number %s',string(iPwd))
        
        case 1
        %utrackPackPaths_cat{iPwd,1} = string(strcat(pwd,filesep,'m',strtrim(movName),filesep,'TrackingPackage'));
        cd(char(utrackPackPaths_cat{iPwd,1})) % safeguarding purpose :: error will be given if path doesn't exist
        utrackPackPaths{iPwd,1} = pwd;
        
        otherwise
        % there should only be 1 folder for each movie; if mergeDir has 
        % multiple cells, it means there are multiple folders that starts
        % with "m-0x-0x" (i.e. "m-01-01" and "m-01-01-beforeDrug"); thus we
        % need to make sure to only go to folder "m-01-01" with uTrack
        % analysis already done on them (that has "TrackingPackage"
        % subfolder.)
        for iPath = 1:length(utrackPackPaths_cat{iPwd,1})
            if isfolder(utrackPackPaths_cat{iPwd,1}(1,iPath))
                cd(utrackPackPaths_cat{iPwd,1}(1,iPath))
                utrackPackPaths{iPwd,1} = pwd;
            end
        end
        
    end % (  switch length(utrackPackPaths_cat{iPwd,1})  )
    
end % (  for iPwd = 1:length(qfsmPackPaths)  )

end
%% ~~~ the end ~~~