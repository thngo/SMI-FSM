function [tassmInl, flagInl, tassmInlFollower1Out, tassmInlFollower2Out] = ...
    retainTassmInliers(tassmMaster, thresholds, tassmFollower1, tassmFollower2)
%RETAINTASSMINLIERS: search for and NaN out outliers in each column of tassm
%
%INPUT:  tassmMaster : (array of n cells), each cell contains 3 cells:
%                    .1st cell contains matrix property for speckle.
%                    .2nd cell contains matrix property for single-molecule.
%                    .(optional) 3rd cell contains matrix for mean position of
%                    single-molecule. 
%
%        thresholds: (1x2 cell) with 1st cell contains vector of threshold
%                    outlier detection for each column in speckle property (i.e.
%                    tassm{i}{1,1}) and 2nd cell contains vector of threshold outlier
%                    detection for each column in single-molecule property (i.e.
%                    tassm{i}{1,2}). If value is NaN, no outlier detection is done on
%                    that column.
%                    Example: thresholds = {[3  NaN  2.5  NaN  4  2.5  2  1.5  2.5] [6  7  2  NaN  NaN  NaN  NaN]};
%
%      tassmFollower1: (optional) similar structure to first input: tassm, but is  
%                    the results from actinPropPerTimeInterval where any
%                    matrix entries with incomplete lifetimes replaced by 
%                    NaN. (see related function cleanActinPropPerTimeInt.m)
%                    Example: tassmClnSpeedLft
%                    If not input, corresponding output is [].
%
%      tassmFollower2: (optional) similar structure to first input: tassm, but is the
%                    results from actinPropPerTimeInterval where any matrix
%                    entries with speeds below threshold replaced by 0.
%                    (see related function cleanActinPropPerTimeInt.m)  
%                    Example: tassmClnSpeed
%                    If not input, corresponding output is [].
%
%OUTPUT: tassmInl  : tassm with inliers only retained. Outliers are NaN out.
%
%        flagInl   : (similar structure with tassm) flag of inliers (1 =
%                    inlier, 0 = outlier)
%
%       tassmInlFollower1Out: similar structure to first output: tassmInl, 
%                    where output flagInl is applied on input tassmClnSpeedLft.
%                    Example: tassmInlClnSpeedLft
%
%       tassmInlFollower2Out: similar structure to first output: tassmInl, 
%                    where output flagInl is applied on input tassmClnSpeed.
%                    Example: tassmInlClnSpeed
%
% Tra Ngo, July 2021
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

%% Input
numCell = length(tassmMaster);
if length(thresholds) ~=2
    error('Input thresholds expects 1x2 cell.');
end

% Check if user only input master tassm (for backward compartibility
% because follower tassm was added later - TN20210716)
if nargin == 2 % (only master tassm was input)
    
    followerTassmCase = 0;
    
elseif nargin == 3 % (if master and only follower 1 were input)
    followerTassmCase = 1;
    % Check that the multiple inputs are of the same tassm size
    numCellFollower1 = length(tassmFollower1);
    if numCellFollower1 ~= numCell
        error('1st and 3rd input must have the same length.')
    end
    
elseif nargin == 4 % (if master and 2 followers were input)
    followerTassmCase = 2;
    % Check that the multiple inputs are of the same tassm size
    numCellFollower1 = length(tassmFollower1);
    numCellFollower2 = length(tassmFollower2);
    if numCellFollower1 ~= numCell || numCellFollower2 ~= numCell
        error('1st, 3rd, and 4th input must have the same length.')
    end
    
else
    
    error('Function expects at least 2 inputs.')
    
end


%% Output:
tassmInl = tassmMaster;
switch followerTassmCase
    case 2 % when 2 followers are input
        tassmInlFollower1Out = tassmFollower1;
        tassmInlFollower2Out = tassmFollower2;
    case 1 % when 1 follower is input
        tassmInlFollower1Out = tassmFollower1;
        tassmInlFollower2Out = [];
    case 0 % when 0 follower is input, only master
        tassmInlFollower1Out = [];
        tassmInlFollower2Out = [];
end

flagInl =  cell(numCell, 1);
for i = 1:numCell
    try
        flagInl{i}{1,1} = true(size(tassmMaster{i}{1,1}));
        flagInl{i}{1,2} = true(size(tassmMaster{i}{1,2}));
    catch
        warning(['Check if tassm of cell ' num2str(i) ' is empty.'])
    end
end

%% Looping through each cell:
for iC = 1:numCell
    
    if isempty(tassmMaster{iC,1}) % skip a cell if its tassm doesn't exist
        continue
    end
    
    % looping through speckle property (iP = 1) & single-molecule property (iP = 2)
    for iP = 1:2
        
        tmpThrLen = length(thresholds{iP});
        tmpMatLen = size(tassmMaster{iC,1}{1,iP},2);
        
        if tmpMatLen ~= tmpThrLen
            error(['Vector length of input thresholds for cell ' num2str(iC) ' needs to be similar to number of column of corresponding tassm matrix.'])
            
        else
            
            % looping through each threshold vector
            for iCurr = 1:tmpThrLen
                
                threshCurr = thresholds{iP}(iCurr);
                
                if isnan(threshCurr) % skip a cell if its threshold is NaN
                    continue
                end
                
                data = tassmMaster{iC,1}{1,iP}(:,iCurr);
                
                % detect outlier
                outlierIdx = isoutlier(data,'quartiles','ThresholdFactor',threshCurr);
                
                % Output inlier flags: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                flagInl{iC}{1,iP}(outlierIdx,iCurr) = false;
                
                % Output tassmInl: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tassmInl{iC}{1,iP}(outlierIdx,iCurr) = NaN;
                
                %% Apply flagInl onto tassmClnSpeedLft, tassmClnSpeed: (NaN out outliers)
                switch followerTassmCase
                    case 2
                        tassmInlFollower1Out{iC}{1,iP}(outlierIdx,iCurr) = NaN;
                        tassmInlFollower2Out{iC}{1,iP}(outlierIdx,iCurr) = NaN;
                    case 1
                        tassmInlFollower1Out{iC}{1,iP}(outlierIdx,iCurr) = NaN;
                end % (if followerTassmFlag)
                
            end % (for iCurr = 1:tmpThrLen)
            
        end % (if tmpMatLen ~= tmpThrLen)
        
    end % (for iP = 1:2)
    
end % (for iC = 1:numCell)

end