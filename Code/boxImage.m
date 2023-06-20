function [boxedImage, boxedImageMeasure] = boxImage(inputImage, boxSize, intstMeasureFlag)
%BOXIMAGE:  split cells into boxes, divide an image up into blocks is by using mat2cell().
%
%SYNOPSIS: [boxedImage, boxedImageIntst] = boxImage(inputImage, boxSize, intstMeasureFlag)
%
%INPUT: inputImage        : (matrix) image with values representing any
%                           measurement (i.e. intensity of speckle, etc.)
%
%       boxSize           : (1x2 vector) [number of pixel row of each box,
%                           number of pixel col of each box] 
%
%       intstMeasureFlag  : (string) indicate whether boxedImageIntst output
%                           'nanmean', 'nansum', or 'sum' of the intensities 
%                           of all the pixels in the box of the inputImage.
%
%OUTPUT: boxedImage       : (cell) each cell lists all measurement within a 
%                           box in the grid.
%
%        boxedImageMeasure: (matrix) sum or mean of measurements of interest.
%
% Previously subfunction within gridIntstCntSpkl.m - TN20221205
% Tra Ngo, Sep 2021
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

% Input control:
if isempty(inputImage)
    error('Input image is empty')
end

% Get the dimensions of the image
[rows, columns, numberOfColorBands] = size(inputImage);

blockSizeR = boxSize(1); % Rows in block.
blockSizeC = boxSize(2); % Columns in block.

% Calculate the size of each block in rows. (blockSizeR & a remainder amount of less than that)
wholeBlockRows = floor(rows / blockSizeR);
if rem(rows, blockSizeR) == 0
    blockVectorR = [blockSizeR * ones(1, wholeBlockRows)];
else
    blockVectorR = [blockSizeR * ones(1, wholeBlockRows), rem(rows, blockSizeR)];
end
% Calculate the size of each block in columns.
wholeBlockCols = floor(columns / blockSizeC);
if rem(columns, blockSizeC) == 0
    blockVectorC = [blockSizeC * ones(1, wholeBlockCols)];
else
    blockVectorC = [blockSizeC * ones(1, wholeBlockCols), rem(columns, blockSizeC)];
end
% divided up the input image into boxes, each pixel info is in a cell in the boxedImage cell array
if numberOfColorBands > 1
    error('Expecting input image to be grayscale!')
    % boxedImage = mat2cell(inputImage, blockVectorR, blockVectorC, numberOfColorBands);
else
    boxedImage = mat2cell(inputImage, blockVectorR, blockVectorC);
end

switch intstMeasureFlag
    case 'nanmean'
        tmp = cellfun(@nanmean, boxedImage,'UniformOutput', false);
        tmp2 = cellfun(@nanmean, tmp,'UniformOutput', false);
        boxedImageMeasure = cell2mat(tmp2);
    case 'sum'
        tmp = cellfun(@sum, boxedImage,'UniformOutput', false);
        tmp2 = cellfun(@sum, tmp,'UniformOutput', false);
        boxedImageMeasure = cell2mat(tmp2);
    case 'nansum'
        tmp = cellfun(@nansum, boxedImage,'UniformOutput', false);
        tmp2 = cellfun(@nansum, tmp,'UniformOutput', false);
        boxedImageMeasure = cell2mat(tmp2);
    otherwise
        error('Boxed image analysis case is not handled.')
end

end