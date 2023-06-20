function [boxedSpeckleList, boxedSpeckleMean, actinPropPerTimeBoxed, spkImg, boxedSpeckleSum, spkIdPerBox] = ...
    gridActinPropPerTimeInt(actinPropPerTimeInt, speckleField, varargin)
% GRIDACTINPROPPERTIMEINT: split movie frame into non-overlapping boxes and get speckle information within each box
%
% SYNOPSIS: [boxedSpeckleList, boxedSpeckleMean, actinPropPerTimeBoxed, spkImg, boxedSpeckleSum, spkIdPerBox] = gridActinPropPerTimeInt(actinPropPerTimeInt, speckleField, varargin)
%
% INPUT:
% Required:
%           actinPropPerTimeInt: (n x 1 cell array of structures) FSM arm
%                             analysis results of SMI-FSM analysis for n
%                             movies; each structure contain speckle
%                             properties for m FSM intervals. For details,
%                             see actinPropPerTimeInterval.m
%
%           speckleField       : (array of char) fieldname from
%                             actinPropPerTimeInt of property of interest.
%            'position'  :  the coordinate is noted by 1 on the output image.
%            'enumerate' :  the coordinate is enumerated by 1, 2, 3, ... in
%                           the listing order in actinStructure.
%                           NOTE: For outputting "spkIdPerBox",
%                           'enumerate' MUST be chosen.
%
% Optional: (as name-value pair)
%           boxWidth  : box width (AKA. number of pixel for each box
%                       row-wise)
%                       Default: 15
%
%           boxHeight : box height (AKA. number of pixel for each box
%                       column-wise)
%                       Default: 15
%
%           nBorder   : number of pixel away from image border that we are
%                       considering the image. This is to correct for dark
%                       pixel around the image border.
%                       Default: 0
%
%           imageSize : image size, determined by the microscope system
%                       (i.e., [341 800] for Olympic Scope, [512 512] for
%                       Nikon Scope, programed on Metomorph).
%                       This can be determined by MovieData.imSize_
%                       Default: [341 800]
%
%           figureFlag: logic flag indicate whether to output figures.
%                       Default: 0 (no figure output)
%
%           maskInfo  : (array of cells) each cell is either a path to
%                       where the mask is stored, or is the mask matrix
%                       itself.
%                       For folder provided, load each mask per frame.
%                       For matrix provided, assume 1 mask for all frames of 1 cell.
%                       Default: [] (no mask, take everything boxes in grid)
%
%
% OUTPUT:   boxedSpeckleList     : (cell array) each cell correspond to a box
%                               in the grid, containing pixels of the grid,
%                               with original image information.
%
%           boxedSpeckleMean     : (matrix) each entry corresponds to a box in
%                               the grid, containing MEAN values of information
%                               from the original image.
%
%           actinPropPerTimeBoxed: (cell array of structure) similar to the
%                               input (actinPropPerTimeInt) but with an
%                               extra field 'gridPosition' indicating the
%                               box index that the speckle fall into within
%                               the grid.
%
%           spkImg               : (cell array of image) each cell is the
%                               original image with speckle information at
%                               speckle position (i.e. speckle intensity
%                               plotted at where speckle is detected,
%                               speckle movement at where speckle is
%                               detected) - independent of grid size.
%
%           boxedSpeckleSum      : (matrix) each entry corresponds to a box in
%                               the grid, containing SUM values of information
%                               from the original image.
%
%           spkIdPerBox          : (cell array) ID of speckle per box,
%                               number of box is dependent on imageSize and
%                               grid size. Each cell contains the following
%                               fields:
%                   .boxId: list the linear index of the boxes, from 1 to
%                           N, N dependent on image size and box size.
%                   .spkId: IDs of speckles that fall within a particular
%                           box. ID is from listing in actinPropPerTimeInt.m
%                   .boxSize: the number of pixel within a particular box.
%                           This is dependent on image size and box size.
%                           Boxes at the image edge will be smaller if
%                           image size is not divisible by box size.
%                   .maskStatus: indicator of whether the box lies
%                           completely within the cell mask (1) or lies
%                           partially at the border of cell mask (2). 
%
% Huong Tra Ngo, Dec 2022
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
ip.FunctionName = 'gridActinPropPerTimeInt';

% Required
ip.addRequired('actinPropPerTimeInt', @iscell) % expect actinPropPerTimeInt to be cell array, each cell correspond to 1 movie
ip.addRequired('speckleField', @ischar) % expect speckleField to be string, matched to fieldnames in actinPropPerTimeInt

% Optional
ip.addParameter('boxWidth', 15, @isscalar)
ip.addParameter('boxHeight', 15, @isscalar)
ip.addParameter('nBorder', 0, @isscalar)
ip.addParameter('imageSize', [341 800], @isvector)
ip.addParameter('figureFlag', 0, @islogical)
ip.addParameter('maskInfo', cell(0,0), @iscell) % mask location or mask; if box is outside of mask => value is NaN
%ip.addParameter('usePartialBox', true, @islogical)

ip.parse(actinPropPerTimeInt, speckleField, varargin{:})

numMov = length(actinPropPerTimeInt);

% (Convenience)
boxWidth = ip.Results.boxWidth;
boxHeight = ip.Results.boxHeight;
nBorder = ip.Results.nBorder;
imageSize = ip.Results.imageSize;
maskInfo = ip.Results.maskInfo; if isempty(maskInfo), maskInfo = cell(numMov,1); end
%usePartialBox = ip.Results.usePartialBox;


mockBox = boxImage(nan(imageSize), [boxWidth boxHeight], 'nanmean'); % create a tmp boxed image
lenBoxList = length(mockBox(:)); clear mockBox

%% Output
boxedSpeckleMean = cell(numMov,0);
boxedSpeckleSum = cell(numMov, 0);
boxedSpeckleList = cell(numMov,0); spkImg = cell(numMov,0);
actinPropPerTimeBoxed = actinPropPerTimeInt;
spkIdPerBox = cell(numMov,0);
numFr = nan(numMov,1);

%% Mask
for kMov = 1:numMov
    numFr(kMov) = size(actinPropPerTimeInt{kMov,1},1);
    for kFr = 1:numFr(kMov)
        
        if ischar(maskInfo{kMov})
            % If mask locations are provided, then go to that folder and load the masks.
            % i.e. maskInfo{1} = /project/biophysics/jaqaman_lab/vegf_tsp1/adasgu/FSM_SMI_ABD_ABDRA/20200113_PD80_LP70_Transfex_point25ug_CD36_C2_D2_day2_halo_3nMJF549_TIMEmngActlowP14/CD36/Actin/m-01-02/QFSMPackage/refined_masks/refined_masks_for_channel_1
            tmp = dir(maskInfo{kMov});
            tmp2 = imread([maskInfo{kMov} filesep tmp(kFr + 2).name]);
            imgMask{kMov,1}{kFr,1} = double(tmp2);
            imgMask{kMov,1}{kFr,1}(imgMask{kMov,1}{kFr,1}==0) = NaN;
        else
            disp('Program assumes user inputs masks for maskInfo');
        end
        
    end % for kFr = 1: length(actinPropPerTimeInt{kMov,1})
end % parfor kMov = 1:numMov

clear kMov kFr

%% Grid each movie
for i = 1:numMov
    actinPropStruct = actinPropPerTimeInt{i,1};
    spkIdPerBox{i,1} = repmat(struct('spkRowId',[],'boxId',[],'spkId',[],'boxSize',[]), numFr(i),1);
    
    % Note: Mask is not needed to be put on the speckle-reconstructed image,
    % because speckle can only come from inside the mask. If we don't
    % care about averaging property by box size, we don't have to worry
    % about boxes being fully inside or only overlapping the mask.
    
    % Step 1: Reconstruct image with information of interest at speckle position
    ogSpkImgs = convSpkCoord2Img(actinPropStruct, imageSize, speckleField);
    
    % spkCntPerBox = c0; spkCntPerBoxFinal = c0; countFinal = c0;
    
    % Step 2: Grid the reconstructed image and calculate measure of
    % interest of speckle counts or or speckle property within each box.
    for iFr = 1:numFr(i) %length(ogSpkImgs)
        
        % Initialize new field:
        spkIdPerBox{i,1}(iFr,1).boxId = (1:lenBoxList)';
        spkIdPerBox{i,1}(iFr,1).spkRowId = cell(lenBoxList,1);
        spkIdPerBox{i,1}(iFr,1).spkId = cell(lenBoxList,1);
        spkIdPerBox{i,1}(iFr,1).boxSize = cell(lenBoxList,1);
        spkIdPerBox{i,1}(iFr,1).maskSize = cell(lenBoxList,1);
        spkIdPerBox{i,1}(iFr,1).maskStatus = cell(lenBoxList,1);
        
        spkImg{i,1}{iFr,1} = ogSpkImgs{iFr}(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
        
        [boxedSpeckleList{i,1}{iFr,1}, boxedSpeckleMean{i,1}{iFr,1}] = ...
            boxImage(spkImg{i,1}{iFr,1}, [boxWidth boxHeight],'nanmean');  % boxedSpeckleMean is for visualization
        % Note that boxedSpeckleList is enumerated by rowId, NOT by actinPropPerTimeInt{i,1}.speckleList.
        
        [~, boxedSpeckleSum{i,1}{iFr,1}] = ...
            boxImage(spkImg{i,1}{iFr,1}, [boxWidth boxHeight],'sum'); % boxedSpeckleSum is for visualization
        
        
        % STEP 2.2: If mask is provided, process it similarly (crop corner + box)
        % and incorporate mask information. (NaN out boxes that are outside of masks)
        
        if ~isempty(imgMask)
            cropImgMask = imgMask{i,1}{iFr,1}(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
            
            [~, partialBoxedMask] = boxImage(cropImgMask, [boxWidth boxHeight],'nanmean'); % Any 0 is box outside of mask.
            [~, boxSize] = boxImage(cropImgMask, [boxWidth boxHeight],'sum');
            [~, boxedMaskArea] = boxImage(cropImgMask, [boxWidth boxHeight],'nansum');
            
            % Incorporate mask into speckle information image
            %             if usePartialBox
            %                 m = partialBoxedMask;
            %             else % do not use partial box
            %                 m = boxSize;
            %                 m(~isnan(m)) = 1;
            %             end
            
            % boxedSpeckleMean and boxedSpeckleSum are for visualization
            boxedSpeckleMean{i,1}{iFr,1} = boxedSpeckleMean{i,1}{iFr,1}.*partialBoxedMask;
            boxedSpeckleSum{i,1}{iFr,1} = boxedSpeckleSum{i,1}{iFr,1}.*partialBoxedMask;
            
            boxedSpeckleList{i,1}{iFr,1}(isnan(partialBoxedMask)) = {NaN};
            % we use partialBox because it keeps the most information of anything related to the
            % mask available. Removing zero matrix and replacing them with NaN is for storage purpose.
            
            % Create map that categorize boxes that are within/partially
            % within/outside mask:
            withinM = double(boxedMaskArea == boxSize); % indicator of boxes within cell mask completely
            tmp = withinM; tmp(:) = 0; tmp(~isnan(partialBoxedMask)) = 1;
            edgeM = tmp ~= withinM; % indicator of boxes within cell mask only partially
            withinM(edgeM == 1) = 2; % adding the cell border indicator (edgeM) to indicators of boxes within cell mask.
            
        end % ( if ~isempty(imgMask) )
        
        % STEP 3: index back into actinPropPerTimeInterval the box that
        % each speckle falls into:
        switch speckleField
            case 'enumerate'
                actinPropPerTimeBoxed{i,1}(iFr).gridPosition = nan(length(actinPropStruct(iFr).speckleList),1);
                for iBox = 1:length(boxedSpeckleList{i,1}{iFr,1}(:))
                    
                    if ~isnan(boxedSpeckleList{i,1}{iFr,1}{iBox}) % Boxes outside the mask are indicated as NaN and are not considered.
                        
                        tmp = find(boxedSpeckleList{i,1}{iFr,1}{iBox}); % Check for existence of speckle for boxes within the mask.
                        
                        if ~isempty(tmp) % Extract speckle information IF THERE EXIST speckle(s) within the checked box.
                            
                            actinPropPerTimeBoxed{i,1}(iFr,1).gridPosition(boxedSpeckleList{i,1}{iFr,1}{iBox}(tmp)) = iBox;
                            %disp(num2str(actinPropPerTimeBoxed{i,1}(iFr,1).gridPosition(boxedSpeckleList{i,1}{iFr,1}{iBox}(tmp))));
                            
                            % Step 4: Update boxes intentified here with corresponding SM identity:
                            spkIdPerBox{i,1}(iFr,1).spkRowId{iBox,1} = boxedSpeckleList{i,1}{iFr,1}{iBox}(tmp);
                            spkIdPerBox{i,1}(iFr,1).spkId{iBox,1} = actinPropPerTimeBoxed{i,1}(iFr,1).speckleList(boxedSpeckleList{i,1}{iFr,1}{iBox}(tmp));
                            
                        end
                        spkIdPerBox{i,1}(iFr,1).boxSize{iBox,1} = boxSize(iBox); %length(boxedSpeckleList{i,1}{iFr,1}{iBox}(:)); % not necessary to recalculate the length!!!
                        spkIdPerBox{i,1}(iFr,1).maskSize{iBox,1} = boxedMaskArea(iBox);
                        spkIdPerBox{i,1}(iFr,1).maskStatus{iBox,1} = withinM(iBox);
                        
                    end
                    
                end
            otherwise
                
        end
        
        
    end % ( for iFr = 1: length(ogSpkImgs) )
    
    
    %% Visualize image
    if ip.Results.figureFlag
        figure, imagesc(boxedSpeckleMean{i,1}{1,1}); title('Speckle properties BOX MEAN, frame 1')
        saveas(gcf,['mean_cell' num2str(i) '.fig'])
        
        figure, imagesc(boxedSpeckleSum{i,1}{1,1}); title('Speckle properties BOX SUM, frame 1')
        saveas(gcf,['sum_cell' num2str(i) '.fig'])
    end
end % ( for i = 1:numMov )
end