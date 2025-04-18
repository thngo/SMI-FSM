function [masksDirOut, maskFolders] = getActinMask(qfsmPackPaths, maskLoc)
%GETACTINMASK: return directory of actual hand-drawn or automated masks that was used for QFSM analysis
%
%INPUT: qfsmPackPaths  : (n x 1 cell) paths leading to QFSM Package results.
%
%       maskLoc        : (n x 1 vector of integer) flag indicating
%                        location of cell ROI mask.
%                            = 0 : importedCellMask; i.e. masks within path 
%                                  provided in qfsmPackPaths (make sure folder 
%                                  name follows a name in extlroi variable.)
%                            = 1 : externalROI or segmentation package's
%                                  mask (for manual user hand-drawn mask).
%                            = 2 : there exists a folder with refined
%                                  masks within the QFSMPackage folder
%                                  provided in qfsmPackPaths (automated
%                                  detection by QFSM package).
%                            = 3 : there exists a folder with refined
%                                  masks within the SegmentationPackage
%                                  folder provided in qfsmPackPaths
%                            If [] or other values, code will crash.
%
%OUTPUT: masksDirOut: (n x 1 cell) masks information. Each cell contains structure
%                     that is the directory of masks.
%        maskFolders: (n x 1 cell) mask folders. This information is
%                     similar to the .folder field of the first output.
%
%Note: this function is used to modularize mask-collecting process in
%actinPropPerTimeInterval.m from SMI-FSM pipeline.
%
%Tra Ngo, Feb 2023
%
% Copyright (C) 2025, Jaqaman Lab - UTSouthwestern 
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
numMov = size(qfsmPackPaths,1);
fs = filesep;
% %naming possible combinations of external masking names which were used in
% %the past. Fourth option may be applied most often as of June 2018
extlroi = ["external_roi","external_mask1","External_mask","External_ROI","external_ROI", ...
    "External_Segmentation"];
% as of February 2019, we have decided that refined_masks from Segmentation
% packages are the same as those from External_ROI, making option 1
% and 3 similar

%% Output
masksDirOut = cell(numMov,0);
maskFolders = cell(numMov,0);

for i = 1:numMov
    if maskLoc(i) == 0
        cd(qfsmPackPaths{i,1})
        cd ../ImportedCellMask
        masks = dir;
        
    elseif ~isempty(qfsmPackPaths{i,1}) && maskLoc(i) ~= 0
        
        %change directory to path to look for folder
        cd(qfsmPackPaths{i,1})
        
        switch maskLoc(i)
            
            case 1
                try
                    % If only External_ROI folder is provided
                    masks = dir([qfsmPackPaths{i,1} fs char(extlroi(isfolder(extlroi)))]);
                catch ME
                    switch ME.identifier
                        case 'MATLAB:catenate:dimensionMismatch'
                            % Windows OS is not case-sensitive and will likely
                            % give multiple possible equivalent folder names.
                            tmpFolder = extlroi(isfolder(extlroi));
                            masks = dir([qfsmPackPaths{i,1} fs char(tmpFolder(1))]);
                        otherwise
                            % These 2 folders contain essentially the same mask
                            % because to use the hand-drawn external mask, we have to pass
                            % that through SegmentationPackage.
                            try
                                masks = dir([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
                            catch
                                % TN20240303: "SegmentationPackage" rebraned as "WindowingPackage"
                                masks = dir([qfsmPackPaths{i,1} fs 'WindowingPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
                            end
                    end
                end
                
            case 2
                masks = dir([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
                
            case 3
                masks = dir([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
                
            otherwise
                disp('error:: Folder of masks does not exist.')
                
        end % switch maskLoc(i)
        
    else
        error('Please input cell mask.')
        
    end % if maskLock(i) == 0 ... elseif ~isempty(qfsmPackPaths{i,1}) && actinConditions.maskLoc(i) ~= 0
    
    masksDirOut{i,1} = masks;
    maskFolders{i,1} = masks(1).folder;

end % for i = 1:numMov

end