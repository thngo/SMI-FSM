function smPropPerTimeIntInMask = retainSmPropInMask(smPropPerTimeIntAll, threshMet, numOutsideMaskFlag)
% RETAINSMPROPINMASK: takes the properties of SMI per FSM interval and the indexes of which SMI tracklet is within mask and discard those tracklets that are outside of the mask.
%
% SYNOPSIS: smPropPerTimeIntInMask = retainSmPropInMask(smPropPerTimeIntAll, threshMet)
%
% EXAMPLE: smPropPerTimeIntInMask = retainSmPropInMask(smPropPerTimeIntSpan2000_tmp, threshMetSpan2000_tmp);
%
% INPUT: smPropPerTimeIntAll: (cell array) SMI properties per interval,
%                             with all tracklets include information both
%                             inside and outside of mask. See
%                             smPropPerTimeInterval.m output for details. 
%
%                  threshMet: (n x 1 cell vectors) contains SM indices that 
%                             exists within cell ROI masks from qFSM.
%
%         numOutsideMaskFlag: (logical) Sanity check flag to display number
%                             of SM being removed due to being outside the
%                             mask. Default: 0 (no display)
%
% OUTPUT: smPropPerTimeIntInMask: similar to input smPropPerTimeIntAll but
%                             containing only SM tracklets within mask.
%
% Tra Ngo, Mar 2023
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

%% Input:
numMov = size(smPropPerTimeIntAll,1);
numMask= size(threshMet,1);
if numMov ~= numMask
    error('Inconsistent inputs. Number of masks and movies must match.')
end
if ~exist('numOutsideMaskFlag','var') || isempty(numOutsideMaskFlag) 
    numOutsideMaskFlag = false;
end


%% Output initialization:
smPropPerTimeIntInMask = cell(numMov,1);

for i = 1:numMov
    
    fnames = fieldnames(smPropPerTimeIntAll{i});
    
    for iFr = 1:size(smPropPerTimeIntAll{i},1) % loop through each interval
        
        for iFl = 1:length(fnames)
            
            switch fnames{iFl}
                case {'smFrames','maskDensity','mergeInfoSpace','splitInfoSpace'}
                    % 'maskDensity' is copied because it is already calculated as number of valid SM tracklets / area of mask
                    % merge and split information is copied because it is independent from SM track as well
                    smPropPerTimeIntInMask{i}(iFr,1).(fnames{iFl}) = smPropPerTimeIntAll{i}(iFr,1).(fnames{iFl});
                otherwise
                    smPropPerTimeIntInMask{i}(iFr,1).(fnames{iFl}) = smPropPerTimeIntAll{i}(iFr,1).(fnames{iFl})(threshMet{i}{iFr},:);
                    
                    % Sanity check:
                    if numOutsideMaskFlag
                        switch fnames{iFl}
                            case{'smList'}
                                disp(['Number of SM outside mask for movie ' num2str(i) ...
                                    ' frame ' num2str(iFr) ': ' num2str(max(smPropPerTimeIntAll{i}(iFr,1).(fnames{iFl})) ...
                                    - length(threshMet{i}{iFr}))])
                        end
                    end

            end
            
        end % for iFl = 1:length(fnames)
        
    end % for iFr = 1:size(smPropPerTimeIntAll{i},1)
    
end % for i = 1:numMov


end