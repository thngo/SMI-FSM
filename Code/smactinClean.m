function [tassmOut, indxConverted] = smactinClean(tassmIn, varName, colIndxArr, thresholdArr)
%SMACTINCLEAN: search for observations in either speckle or single-molecule properties in tassm that are LOWER OR EQUAL TO input threshold and replace them with NaN.
%
% SYNOPSIS: [tassmOut, indxConverted] = smactinClean(tassmIn, varName, colIndx, threshold)
%
% Example: [tassmOut, indxConverted] = smactinClean(tassmIn, 'speckle', 1, 0.01)
%
% INPUT:    tassmIn: (array of cells) each cell contains 3 cells denoting
%                    matrix for properties of speckle, SM, SM position.)
%
%           varName: (string) either exactly as:
%                       'speckle': function looks into first cell of
%                       tassm{i}, containing speckle properties.
%                       'single molecule': function looks into 2nd cell of
%                       tassm{i}, containing SM properties.
%
%           colIndxArr: (vector of integer) the colume index of the
%                       property of interest 
%                       Example: colIndxArr = [1 9]; (for speckle speed and 
%                                speckle average speed) 
%
%           thresholdArr: (vector of double aka. real number) the threshold 
%                       where value equal or below which is converted to NaN.
%                       Example: thresholdArr = [0 0];
%
% OUTPUT:   tassmOut: (array of cells) each cell contains 3 cells denoting
%                     matrix for properties of speckle, SM, SM position.)
%
%           indxConverted: (array of cells that contain n-by-1 vector of
%                          integer) 1 indicates that that observation was
%                          below threshold and converted to NaN in
%                          tassmOut.
%
% Tra Ngo, Jun 2021
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

%% Input control & variable initiations
numCell = length(tassmIn);
numProp = length(colIndxArr);

switch varName
    case 'speckle'
        iProp = 1;
    case 'single molecule'
        iProp = 2;
    otherwise
        error('Unexpected input. "varName" expects only either "speckle" or "single molecule" string.')
end

% Control: length of colIndxArr has to equal length of thresholdArr
if length(thresholdArr) ~= numProp
    error('colIndxArr and thresholdArr inputs have to be vectors of the same length.')
end

%% Output
tassmOut = tassmIn;
indxConverted = cell(numCell,1);

for iCol = 1:numProp % loop through each properties of interest
    
    % grab colume index and according threshold to process
    colIndx = colIndxArr(iCol);
    threshold = thresholdArr(iCol);
    
    for iCell = 1:numCell % loop through each cell in tassmIn
        
        % skip empty cell
        if isempty(tassmIn{iCell,1})
            continue
        end
        
        % grab colume containing all observations for property of interest from
        % current cell:
        allObs = tassmIn{iCell,1}{1,iProp}(:,colIndx);
        
        % find index of observations that does not meet the input threshold:
        indxNan = (allObs <= threshold);
        
        % replace the observations that does not meet the threshold with NaN:
        tassmOut{iCell,1}{1,iProp}(indxNan,colIndx) = NaN;
        
        % output index of observatioins that does not meet input threshold for
        % current cell:
        indxConverted{iCell,1}{1,iCol} = indxNan;
        
    end
end

end