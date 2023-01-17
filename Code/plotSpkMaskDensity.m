function [maskDens, maskDensMat, maskDensMatCov] = plotSpkMaskDensity(actinPropPerTimeStructArr, subgroupArr, pixelSize)
%PLOTSPKMASKDENSITY collect and plot speckle Mask Density
%
%SYNOPSIS: [maskDens, maskDensMat, maskDensMatCov] = plotSpkMaskDensity(actinPropPerTimeStructArr, subgroupArr, pixelSize)
%
%INPUT: actinPropPerTimeStructArr: (cell array of structure) output of
%                                  actinPropPerTimeInt.m
%
%       subgroupArr              : (vector) each entry is index number of
%                                  cell to analyse from actinPropPerTimeStructArr
%
%       pixelSize                : (1x2 cell array, size of each pixel and its unit. 
%                                  . empty : the output density do not get converted
%                                           (unit = speckle per pixel^2)
%                                  . [any number, unit string] : the
%                                           output density get converted to
%                                           unit indicated in unit string.
%                                       i.e.: {90, 'nm'}
%
%OUTPUT: maskDens                 : (cell array) each cell is a cell from actinPropPerTimeStructArr
%
%        maskDensMat              : (matrix) row = frame, col = cell
%
%        maskDensMatCov           : (matrix) similar to maskDensMat but in
%                                   converted unit.
%
%Tra Ngo, Nov 2021
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

maskDensMatCov = [];
emptyCount = 1;

for i = 1:length(subgroupArr)
    try
    maskDens{i,1} = [actinPropPerTimeStructArr{subgroupArr(i),1}(1:end).speckleMaskDensity]';
    catch
        warning('If input contains an empty cell, that cell will be removed from output matrix!')
        if isempty(actinPropPerTimeStructArr{subgroupArr(i),1})
            emptyIndex(emptyCount) = i;
            emptyCount = 1 + emptyCount;
        end
    end
end

maskDensMat = cell2mat(maskDens'); % row = frame; col = cell

if exist('emptyIndex','var') && length(emptyIndex) == 1 % if input contain an empty cell, add NaN so that the indexing is maintained.
    nRow = size(maskDensMat,1);
    maskDensMat = [maskDensMat(:,1:emptyIndex-1) nan(nRow,1) maskDensMat(:,emptyIndex:end)];
elseif exist('emptyIndex','var') && length(emptyIndex) ~=1 
    error('Code note currently written to handle multiple empty cells');
end

if isempty(pixelSize)
    unitStr = 'pixel'; % Get the unit the density will be in
    
    figure, plot(maskDensMat,'DisplayName','maskDensMat')
else
    unitStr = pixelSize{2}; % Get the unit the density will be in
    area1px = pixelSize{1}^2;
    maskDensMatCov = maskDensMat/area1px; % (converted density) = (old density) / (new equivalent area)
    
    figure, plot(maskDensMatCov,'DisplayName','maskDensMat')
end

xlabel('frames'); ylabel(['mask density (speckle per ' unitStr '^2)']); title('Speckle mask density at each frame for each cell')


end

