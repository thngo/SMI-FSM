function [actinPropPerTimeIntClnSpeed,actinPropPerTimeIntClnLft,...
    actinPropPerTimeIntClnSpeedLft] = cleanActinPropPerTimeInt(...
    actinPropPerTimeInt,speedThresh)
%CLEANACTINPROPPERTIMEINT replaces small speckle speeds and incomplete speckle lifetimes with 0 and NaN, respectively
%
%SYNOPSIS [actinPropPerTimeIntClnSpeed,actinPropPerTimeIntClnLft,...
%     actinPropPerTimeIntClnSpeedLft] = cleanActinPropPerTimeInt(...
%     actinPropPerTimeInt,speedThresh)
%
%INPUT  actinPropPerTimeInt: Cell array containing structures with various
%                            fields containing actin speckle properties, as
%                            output by actinPropPerTimeInterval.
%       speedThresh        : Threshold below which speeds are replaced with
%                            NaN (because they are difficult to tell apart
%                            from apparent movement due to noise).
%                            Example: speedThresh = 1.5;
%
%OUTPUT actinPropPerTimeIntClnSpeed: Same as input actinPropPerTimeInt, but
%                            with speeds below threshold replaced by 0.
%       actinPropPerTimeIntClnLft: Same as input actinPropPerTimeInt, but
%                            with incomplete lifetimes replaced by NaN.
%       actinPropPerTimeIntClnSpeedLft: Same as input actinPropPerTimeInt, but
%                            with speeds below threshold replaced by 0 and
%                            incomplete lifetimes replaced by NaN.
%
%Khuloud Jaqaman, July 2021
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

%initialize output
[actinPropPerTimeIntClnSpeed,actinPropPerTimeIntClnLft,...
    actinPropPerTimeIntClnSpeedLft] = deal(actinPropPerTimeInt);

%loop over movies
numMovies = length(actinPropPerTimeInt);
for iMovie = 1 : numMovies
    
    %get current movie info from input cell array
    actinPropCurr = actinPropPerTimeInt{iMovie};
    [actinPropClnSpeed,actinPropClnLft,actinPropClnSpeedLft] = deal(actinPropCurr);
    
    %loop over frames
    numFrames = length(actinPropCurr);
    for iFrame = 1 : numFrames
        
        %speckles with speed below threshold
        speckleSpeed = actinPropCurr(iFrame).speckleSpeed;
        indxBad = find(speckleSpeed < speedThresh);
        
        actinPropClnSpeed(iFrame).speckleSpeed(indxBad) = 0;
        actinPropClnSpeed(iFrame).speckleVelocity(indxBad,:) = 0;
        actinPropClnSpeed(iFrame).speckleMvmtCohere(indxBad) = NaN;
       
        actinPropClnSpeedLft(iFrame).speckleSpeed(indxBad) = 0;
        actinPropClnSpeedLft(iFrame).speckleVelocity(indxBad,:) = 0;
        actinPropClnSpeedLft(iFrame).speckleMvmtCohere(indxBad) = NaN;
        
        %speckles with incomplete lifetime
        speckleLftFlag = actinPropCurr(iFrame).flagCompleteLft;
        indxBad = find(speckleLftFlag==0);
        
        actinPropClnLft(iFrame).speckleLifetime(indxBad) = NaN;
        
        actinPropClnSpeedLft(iFrame).speckleLifetime(indxBad) = NaN;
        
    end %(for iFrame = 1 : numFrames)
    
    %store in output cell array
    actinPropPerTimeIntClnSpeed{iMovie} = actinPropClnSpeed;
    actinPropPerTimeIntClnLft{iMovie} = actinPropClnLft;
    actinPropPerTimeIntClnSpeedLft{iMovie} = actinPropClnSpeedLft;
    
end %(for iMovie = 1 : numMovies)

end

%% ~~~ the end ~~~