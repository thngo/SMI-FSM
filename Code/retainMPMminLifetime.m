function [mpmMinLft, lftimeMat, speckleLftTot]= retainMPMminLifetime(mpm,minLft)
%RETAINMPMMINLIFETIME calculate lifetime and replace speckle position of speckles whose lifetime is less than a given number with zeros in magic position matrix MPM
%
%SYNOPSIS mpmMinLft = retainMPMminLifetime(mpm,minLft)
%
%INPUT      mpm    : speckle track matrix (Magic Position Matrix) as output
%                    by qFSM software
%                    e.g. load('.../QFSMPackage/speckleTracks/_tracks.mat','MPM');
%
%           minLft : (integer) minimum lifespan of the speckle (i.e. minimum
%                   number of frames a speckle has to exist) to retain in MPM
%                           
%OUTPUT  mpmMinLft :    matrix that replicates Magic Position Matrix MPM with only  
%                       speckles that exist at least a "minLft" number of frame
%
%        lftimeMat :    matrix containing lifetime of each speckle
%                       corresponding to the MPM
%
%        speckleLftTot: (BUGGY, DUPLICATE SPECKLE AT EACH FRAME, DO NOT
%                       USE) horizontal vector that contains lifetime of
%                       all speckle in input MPM.
% 
%Tra H. Ngo, February 2019; added lftimeMat Aug. 2019
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

mpmMinLft = mpm;
speckleLftTot = [];
% Extract only the x-position of speckle out of MPM
mpmX = mpm(:,1:2:end);

numRow = size(mpmX,1);

% Create matrix mpmXL of double with 1 indicating existence of speckle and 
% 0 indicating death of speckle
mpmXL = [mpmX>0 zeros(numRow,1)];
% create lifetime matrix correspond to MPM
lftimeMat = zeros(size(mpm));
for iRow = 1 : numRow
    
    iCol = find(mpmXL(iRow,:)==0); % extract indexes of zero entries row-wise in logic matrix mpmXL
    speckleLft = diff([0 iCol])-1; % speckle lifetime row-wise of speckles in MPM
    
    speckleLftTot = horzcat(speckleLftTot, speckleLft); % concatenate calculated speckle lifetime from each row
    
    
    for jCol = 1:length(iCol)
        % Replace position of speckle whose lifetime is less than minLft number
        % with zeros in mpmX matrix
        if speckleLft(jCol) < minLft
            for subtractor = 1:(2*speckleLft(jCol))
                mpmMinLft(iRow, 2*iCol(jCol) - 1 - subtractor) = 0;
            end
        end
        % fill-in lifetime matrix correspond to MPM
        for subtractor = 1:(2*speckleLft(jCol))
            lftimeMat(iRow, 2*iCol(jCol) - 1 - subtractor) = speckleLft(jCol);
        end
    end   
    
    %KJ: This would be an alternative that minimizes the use of for loops and if statements
    %    I thought this would be faster, but apparently it is not
    %    Maybe it will be faster for much bigger matrices
    %
    %     %find speckles with lifetime less than minLft
    %     iSpeckleShort = find(speckleLft > 0 & speckleLft < minLft);
    %
    %     % Replace position of speckle whose lifetime is less than minLft number
    %     % with zeros in mpmX matrix
    %     for jCol = iSpeckleShort
    %         subtractor = 1:(2*speckleLft(jCol));
    %         mpmMinLft(iRow,2*iCol(jCol) - 1 - subtractor) = 0;
    %     end
    
end