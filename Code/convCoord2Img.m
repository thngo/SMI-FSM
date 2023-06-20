function [outputImg,outputImgArray] = convCoord2Img(coord, imageSize, imageType, bgType)
% CONVCORD2IMG convert x-,y- positions to a 2D image
%
% SYNOPSIS: outputImages = convCoord2Img(coord, imageSize, imageType, bgType)
%
% INPUT: coord     : (matrix) x-,y-positions, with x in 1st col, y in 2nd col
%
%        imageSize : (1x2 double) input that provide information on image size.
%
%        imageType : (string)
%            'position'  :  the coordinate is noted by 1 on the output image.
%            'enumerate' :  the coordinate is enumerated by 1, 2, 3, ... in
%                           the listing order in coord.
%
%         bgType   : (string) background type
%            'nan' : the empty space is denoted by NaNs
%            'zeros' : the empty space is denoted by 0s
%
% OUTPUT: outputImg: (matrix of size imageSize) 0s denoting blank space,
%                    non-zero numbers denothing input coordinations.
%
%         outputImgArray: (cell) created to handle sub-pixel localization
%                    of SM.
%
% Written based on convSpkCoord2Img.m
% Tra Ngo, Feb 2023
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

%% Output
switch bgType
    case 'nan'
        outputImg = nan(imageSize);
    case 'zeros'
        outputImg = zeros(imageSize); % generate a blank image
    otherwise
        error('Invalid background type. Expect nans or zero')
end
outputImgArray = cell(0,0);

%% Get boxed image:
indexTest = sub2ind(imageSize, coord(:,2), coord(:,1));

switch imageType
    case 'position'
        outputImg(indexTest) = 1; % turn all positions of speckle into 1
    case 'enumerate'
        iCnt = 1;
        for iEnumerate = 1:length(indexTest) % enumerate all positions into 1:1:length(indexTest)
            if isnan(outputImg(indexTest(iEnumerate))) || outputImg(indexTest(iEnumerate)) == 0
                outputImg(indexTest(iEnumerate)) = iEnumerate ;
            else
                outputImgArray{iCnt} = outputImg;
                outputImgArray{iCnt}(indexTest(iEnumerate)) = iEnumerate ;
                iCnt = iCnt + 1;
            end
        end
    otherwise
        error('Field not handled correctly (wrong size or wrong field name).')

end % ( switch imageType)


end