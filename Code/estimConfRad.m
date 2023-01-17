function [confRadTmp,centerTmp] = estimConfRad(tracks,probDim,confRadMin)
%ESTIMCONFRAD estimate confinement radius of tracklets
%
%INPUT:         tracks (.tracksCoordAmpCG)
%                           Matrix indicating the positions and amplitudes
%                           of the tracked features to be plotted. Number
%                           of rows = number of tracks, while number of
%                           columns = 8*number of time points. Each row
%                           consists of
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           NaN is used to indicate time points where a
%                           track does not exist.
%                      
%
%               probDim: Problem dimensionality.
%
%               confRadMin: if = 1, min(eigenVal) is used
%                           if = 0, mean(eigenVal) is used (default = 0 in basicTransientDiffusionAnalysisv1.m)
%
%OUTPUT: confRadTmp: Confinement/Diffusion Radius
%
%        centerTmp:  Center of confinement/diffusion neighborhood
%
%From basicTransientDiffusionAnalysisv1.m subfunction by Tony Vega July 2016
%Moved Dec 2019 Tra Ngo
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

    %get subpart's coordinates
    xCoord = tracks(1:8:end)';
    yCoord = tracks(2:8:end)';
    zCoord = tracks(3:8:end)';
    xyzCoord = [xCoord yCoord zCoord];
    
    if all(isnan(xyzCoord))
        
        %handle case of all NaNs
        confRadTmp = NaN;
        centerTmp = NaN;
        
    else
        
        %find the eigenvalues and eigenvectors of the variance-covariance
        %matrix of this track's positions
        eigenVal = eig(nancov(xyzCoord(:,1:probDim)));
        
        %calculate the track's confinement radius
        if confRadMin
            confRadTmp = sqrt( min(eigenVal) * (probDim + 2) );
        else
            confRadTmp = sqrt( mean(eigenVal) * (probDim + 2) );
        end
        
        %calculate the track's center
        centerTmp = nanmean(xyzCoord(:,1:probDim));
        
    end
    
end