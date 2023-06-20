function [edgeClass, gridClassInfo, imgErodedMask] = classifyMaskGridByErosion(mask, gridSize, numErosion, se, nBorder, senseCheck)
%CLASSIFYMASKGRIDBYEROSION: grid mask into boxes and classify those boxes into edge class and center class based on iterative erosion
%
%SYNOPSIS: [edgeClass, gridClassInfo, imgErodedMask] = classifyMaskGridByErosion(mask, gridSize, numErosion, se, nBorder, senseCheck)
%
%INPUT: mask      : binary image of cell mask
%
%       gridSize  : (1x2 vector matrix) width and height of each box that will
%                   grid out the cell mask image.
%
%       numErosion: (integer) Number of erosion that will erode the
%                   cell mask edge.
%
%       se        : structuring elements for mask erosion
%
%       nBorder   : (integer) number of pixel around mask border to remove.
%
%       senseCheck: (logical) Flag to plot the cell mask at each
%                   erosions.
%
%OUTPUT: edgeClass: (n*m x 1 vector matrix)
%                   0: correspond to background
%                   1: all boxes inside the mask of the last eroding-iteration as "center"
%                   2: all boxes that have been edge box (through all iteration) as "edge"
%
%        gridClassInfo: contains all information for classification of each
%                       box in the grid, at each erosion iteration.
%                       Contains the folling fields:
%                       .boxSize: number of pixels within a box.
%                       .maskSize: number of pixels within a box that is
%                                   also within the cell mask.
%                       .maskStatus
%                       .maskStatusAllErode: the erosion iteration at
%                                   which this box is classified as "edge".
%                            2 = 1st cell edge (from original image)
%                            3 = 2nd cell edge (from 1st erosion)
%                            4 = 3rd cell edge (from 2st erosion)
%
%       imgErodedMask: (n x 1 cell array, n = number of erosion) cell mask
%                      after each erosion iteration.
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

initCell = cell(1,(numErosion+1));
partialBoxedMask = initCell; boxSize = initCell; boxedMaskArea = initCell;
withinM = initCell; edgeM = initCell; withinM_tmp = initCell;

if ~isempty(mask)
    
    imgErodedMask = cell(numErosion,1);
    imgErodedMaskNumeric = cell(numErosion,1);
    
    imgMask = double(mask);
    imgMask(imgMask==0) = NaN;
    
    imgErodedMask{1} = mask; % iteration 1 (before eroding) = real mask
    imgErodedMaskNumeric{1} = imgMask; % iteration 1 (before eroding) = real mask
    
    % Iterative erosion
    for m = 1:(numErosion)
        imgErodedMask{m+1} = imerode(imgErodedMask{m}, se); % store new eroded mask
        imgErodedMaskNumeric{m+1} = double(imgErodedMask{m+1}); % convert eroded mask to double
        imgErodedMaskNumeric{m+1}(imgErodedMaskNumeric{m+1} == 0) = NaN; % convert 0 to NaN in eroded mask
    end
    
    % Calculate mask size and box size and mask status
    for j = 1:(numErosion+1)
        cropImgMaskTmp = imgErodedMaskNumeric{j}(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
        
        % Since we do not change the grid, we only need to calculate size of each box once at the start. 
        if j == 1 
            [~, boxSize] = boxImage(cropImgMaskTmp, gridSize,'sum');
        end
        
        [~, partialBoxedMask{j}] = boxImage(cropImgMaskTmp, gridSize,'nanmean'); % Any NaN (or 0) is box outside of mask.
        [~, boxedMaskArea{j}] = boxImage(cropImgMaskTmp, gridSize,'nansum'); % Any 0 is box outside of mask.
        
        % Create map that categorize boxes that are within/partially within/outside mask:
        withinM{j} = double(boxedMaskArea{j} == boxSize); % indicator of boxes within cell mask completely
        tmp = withinM{j}; tmp(:) = 0; tmp(~isnan(partialBoxedMask{j})) = 1; % indicator of boxes that touch (fully + partially) any part of the cell mask
        edgeM{j} = tmp ~= withinM{j}; % indicator of boxes within cell mask only partially
        withinM{j}(edgeM{j} == 1) = j+1; % adding the cell border indicator (edgeM) to indicators of boxes within cell mask.
        
        withinM_tmp{j} = withinM{j}(:);
    end % for j = 1:(numErosion+1)
    
    withinM_tmp_all = max(cell2mat(withinM_tmp),[],2); % BECAREFUL WHEN MAX IS 1, IS IT ALWAYS 1 OR IS IT [1 0 0 0 0]
    
    for i = 1:(numErosion+1)
        if i == 1
            gridClassInfo.boxSize = boxSize;
            gridClassInfo.maskSize = boxedMaskArea{i};
            gridClassInfo.maskStatus = withinM{i};
            gridClassInfo.maskStatusAllErode = withinM_tmp_all;
        else
            gridClassInfo.(['maskSize' num2str(i)]) = boxedMaskArea{i};
            gridClassInfo.(['maskStatus' num2str(i)]) = withinM{i};
        end % if i == 1
    end % for i = 1:(numErosion+1)
    
    if senseCheck
        
        [row, col] = ind2sub(size(gridClassInfo.boxSize), 1:length(gridClassInfo.maskStatusAllErode));
        
        indx = gridClassInfo.maskStatusAllErode ~= 0;
        for i = 1:length(gridClassInfo.maskStatusAllErode)
            statusTmp{i} = num2str(gridClassInfo.maskStatusAllErode(i));
        end
        
        figure;
        nCol = ceil(sqrt(numErosion+1));
        nRow = nCol - (nCol^2 - (numErosion + 1) > nCol - 1);
        for k = 1:(numErosion+1)
            subplot(nRow, nCol, k)
            imshow(imgErodedMask{k}); hold on;
            text(col(indx)*(size(mask,2)/size(gridClassInfo.boxSize,2))-(gridSize(2)/1.5),... % decrease col => text shift left
                row(indx)*(size(mask,1)/size(gridClassInfo.boxSize,1))+(gridSize(1)/8),... % increase row => text shift down
                statusTmp(indx),'Color', 'r', 'FontSize',7);
            pause(3)
        end
        
    end % if senseCheck
    
end % if ~isempty(mask)

edgeClass = zeros(numel(gridClassInfo.boxSize),1); % initialize all boxes with 0, correspond to background
edgeClass(withinM_tmp_all == 1) = 1; % denote all boxes inside the mask of the last eroding-iteration as "center" (1)
edgeClass(withinM_tmp_all > 1) = 2; % denote all boxes that have been edge box (through all iteration) as "edge" (2)

end