function [maskListClass,imgErodedMask] = classifyMaskGrid(maskFolderList, varargin)
%CLASSIFYMASKGRID: takes cell masks and return classifications of cell edge and center regions per box in user-defined grid
%
%SYNOPSIS: [maskListClass,imgErodedMask] = classifyMaskGrid(maskFolderList, varargin)
%
%INPUT: Required. maskFolderList: list of folder with mask files as .tif files
%
%       Optional.
%           numErosion: (integer) Number of erosion that will erode the
%                       cell mask edge. Default = 0.
%
%           se        : structuring elements for mask erosion. Default =
%                       strel('rectangle', [15 15]).
%
%           nBorder   : (integer) number of pixel around mask border to
%                       remove. Default = 0.
%
%       	boxWidthHeight: (vector) width and height of each box that will
%                       grid out the cell mask image. Default = [15 15].
%
%           senseCheck: (logical) Flag to plot the cell mask at each
%                       erosions. Default = false. 
%
%OUTPUT: maskListClass: (cell array) information for each box of the grid,
%                       contain the following fields:
%                       .edgeClass = flag indicating whether a
%                                   particular box correspond to background
%                                   (0), center subcellular region (1), or
%                                   edge subcellular region (2).
%                       .maskStatusAllErode = the erosion iteration at
%                                   which this box is classified as "edge".
%                       .boxSize = number of pixels within a box.
%                       .maskSize = number of pixels within a box that is
%                                   also within the cell mask.
%                       .maskStatus(iBox,1) = gridClassInfo.maskStatus(iBox);
%        imgErodedMask: (cell array) cell border indicator (edgeM) to
%                       indicators of boxes within cell mask totally , vs
%                       within cell mask partially
%
%Tra Ngo, April 2023
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
% along with SMI-FSM. If not, see <http://www.gnu.org/licenses/&gt;.
% 

%% Input
ip = inputParser;
ip.CaseSensitive = false;
ip.FunctionName = 'classifyMaskGrid';

% Required
ip.addRequired('maskFolderList', @iscell) % expect cell array of directories

% Optional
ip.addParameter('numErosion', 0, @isnumeric)
ip.addParameter('se', strel('rectangle', [15 15]), @isscalar)
ip.addParameter('nBorder', 0, @isscalar)
ip.addParameter('boxWidthHeight', [15 15], @isvector)
ip.addParameter('senseCheck', false, @islogical) % false => no output image

ip.parse(maskFolderList, varargin{:})

% (Convenience)
numErosion = ip.Results.numErosion;
se = ip.Results.se;
nBorder = ip.Results.nBorder;
boxWidthHeight = ip.Results.boxWidthHeight;
senseCheck = ip.Results.senseCheck;


%% Mask, input, variable initialization:
disp(['Classifying grid of ' num2str(boxWidthHeight(1)) 'x' num2str(boxWidthHeight(2)) ' boxes with ' num2str(numErosion) ' erosions, ignoring ' num2str(nBorder) ' px at the image border.']);

numMov = length(maskFolderList);
numFr = nan(numMov,1);
imgMask = cell(numMov,0);
lenBoxList = cell(numMov,1);
maskListClass = cell(numMov,1);
imgErodedMask = cell(numMov,1);

for iMov = 1:numMov
    
    tmp = dir(maskFolderList{iMov});
    numFr(iMov) = size(tmp,1)-2;
    
    for kFr = 1:numFr(iMov)
        
        if ischar(maskFolderList{iMov})
            % If mask locations are provided, then go to that folder and load the masks.
            loadedMask{iMov,1}{kFr,1} = imread([maskFolderList{iMov} filesep tmp(kFr + 2).name]);
            imgMask{iMov,1}{kFr,1} = double(loadedMask{iMov,1}{kFr,1});
            imgMask{iMov,1}{kFr,1}(imgMask{iMov,1}{kFr,1}==0) = NaN;
            
            mockBox = boxImage(nan(size(loadedMask{iMov,1}{kFr,1})), boxWidthHeight, 'nanmean'); % create a tmp boxed image
            lenBoxList{iMov,1}(kFr,1) = length(mockBox(:)); clear mockBox
            
        else
            error('Input maskFolderList needs to have characters in each entry');
        end
        
    end % for kFr = 1: length(actinPropPerTimeInt{iMov,1})
end % parfor iMov = 1:numMov

clear iMov kFr

%% Grid each movie
for i = 1:numMov
    
    maskListClass{i,1} = repmat(struct('boxId',[],'edgeClass',[],'maskStatusAllErode',[]), numFr(i),1);
    
    for iFr = 1:numFr(i)  % looping through each FSM frame
        
        % Initialize new field:
        maskListClass{i,1}(iFr,1).boxId = [1:lenBoxList{i,1}(iFr,1)]';
        
        % STEP 0: Process mask information
        if (iFr == 1) || (iFr == numFr(i))
            [edgeClass, gridClassInfo, imgErodedMask{i,1}{iFr,1}] = classifyMaskGridByErosion(loadedMask{i,1}{iFr,1}, boxWidthHeight, numErosion, se, nBorder, senseCheck);
        else
            [edgeClass, gridClassInfo, imgErodedMask{i,1}{iFr,1}] = classifyMaskGridByErosion(loadedMask{i,1}{iFr,1}, boxWidthHeight, numErosion, se, nBorder, false);
        end
        
        
        % initialize the fields:
        maskListClass{i,1}(iFr,1).boxSize    = nan(lenBoxList{i,1}(iFr,1),1);
        maskListClass{i,1}(iFr,1).maskSize   = nan(lenBoxList{i,1}(iFr,1),1);
        maskListClass{i,1}(iFr,1).maskStatus = nan(lenBoxList{i,1}(iFr,1),1);
        
        for iBox = 1:lenBoxList{i,1}(iFr,1)
            maskListClass{i,1}(iFr,1).edgeClass = edgeClass;
            maskListClass{i,1}(iFr,1).maskStatusAllErode = gridClassInfo.maskStatusAllErode;
            maskListClass{i,1}(iFr,1).boxSize(iBox,1) = gridClassInfo.boxSize(iBox);
            maskListClass{i,1}(iFr,1).maskSize(iBox,1) = gridClassInfo.maskSize(iBox);
            maskListClass{i,1}(iFr,1).maskStatus(iBox,1) = gridClassInfo.maskStatus(iBox);
        end % ( for iBox = 1:lenBoxList{i,1}(iFr,1) )
        
    end % for iFr = 1:numFr  % looping through each FSM frame
    
    
end

end