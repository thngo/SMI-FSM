function smPropPerTimeIntByBoxClass = stratifySmTrackletByBoxClass(smPropPerTimeInterval, boxClass, boxClassCondition, smTrackletIdPerBox)
%STRATIFYSMTRACKLETBYBOXCLASS: split the SM data into groups based on their box classification
%
%SYNOPSIS: smPropPerTimeIntByBoxClass = stratifySmTrackletByBoxClass(smPropPerTimeInterval, boxClass, boxClassCondition, smTrackletIdPerBox)
%
%INPUT: smPropPerTimeInterval: (n x 1 cell array of structures) SMI arm
%                            analysis results of SMI-FSM analysis for n
%                            movies; contains various properties of SM; see
%                            smPropPerTimeInterval.m and gridSmPerTimeInt.m
%                            for details. 
%
%       boxClass: (cell array) information for each box of the grid,
%                            including classification information. See
%                            output of classifyMaskGrid.m for details.  
%
%       boxClassCondition: (structure) with following fields:
%           .extractClass: (1D vector) integer encoding of the box classes
%                         user wants to extract. For example:
%                         If input [1], user wants to extract only boxes
%                         classified as "1".
%                         If input [1 2], user wants to extract boxes
%                         classified as "1" and boxes classified as "2",
%                         which will be saved as 2 separate cells in the
%                         output array.
%           .firstFrame : the first frame of mask that should be considered
%                         box classification of SM movie. (often similar to
%                         .firstFsmFrame of smConditions, i.e. = 2 (for
%                         darkframe-scheme dataset)).
%
%       smTrackletIdPerBox: (cell array of structure) ID of SM per box,
%                          number of box is dependent on imageSize and grid
%                          size. See output of gridSmPerTimeInt.m for details.
%
%OUTPUT: smPropPerTimeIntByBoxClass: (n x m cell array) contains various
%                         properties of SM, splitted into groups based on
%                         box classification. 
%                         n = number of movie from input smPropPerTimeInterval; 
%                         m = number of box classes user wants to extract.
%
% Tra Ngo April 2023
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

numMov = size(smPropPerTimeInterval,1);
numMask = size(boxClass,1);
numBoxId = size(smTrackletIdPerBox,1);
if numMov ~= numMask || numMov ~= numBoxId
    error('The movie array inputs must have similar length.')
end

numClass = length(boxClassCondition.extractClass);
smPropPerTimeIntByBoxClass = cell(numMov,numClass);

if max(boxClassCondition.extractClass) > 2
    warning('Function uses assumption of class 1 (center) & 2 (edge) for classification.')
end

for i = 1:numMov
    numFr = size(smPropPerTimeInterval{i},1);
    
    for iFr = 1:numFr
        
        iMaskFr = iFr + boxClassCondition.firstFrame - 1;
        
        tmp_sm = rmfield(smPropPerTimeInterval{i}(iFr), {'smFrames', 'maskDensity'});
        fieldList = fieldnames(tmp_sm);
        
        tmp_numTrk = length(tmp_sm.boxAssigned);
        % STEP 1: For current interval, get index of each class. (Read the
        % vector indicating which box localization of SM is in and assign
        % the types to the box. Final output is the type that the SM
        % tracklet sits in majority of the time.) -------------------------
        tmp_smIndxClass = false(length(tmp_sm.boxAssigned), numClass);
        tmp_percInClassAssigned = zeros(length(tmp_sm.boxAssigned), numClass);
        for iTrk = 1:tmp_numTrk
            
            tmp_class = boxClass{i}(iMaskFr).edgeClass(tmp_sm.gridPosition{iTrk, 1});
            
            % Get the 1 box type that track this matched to most of the time
            uniqueClass = unique(tmp_class);
            numLocalization = length(tmp_class);
            percInBox = nan(length(uniqueClass),1);
            for iUniq = 1: length(uniqueClass) % loop through each class and count the percentage of being in each class.
                percInBox(iUniq) = nansum(tmp_class == uniqueClass(iUniq))/numLocalization;
            end
            indxBoxAssigned = percInBox == max(percInBox); % index of the box that SM is being assigned to.
            % ASSUMPTION: class 2 is dominant. If SM has both class 1 &
            % 2, then its final assignment is class 2.
            if sum(indxBoxAssigned) == 1
                classAssigned = uniqueClass(indxBoxAssigned);
                tmp_percInClassAssigned(iTrk,boxClassCondition.extractClass == classAssigned) = percInBox(indxBoxAssigned);
            elseif sum(indxBoxAssigned) > 1
                classAssigned = 2;
                tmp_percInClassAssigned(iTrk,boxClassCondition.extractClass == classAssigned) = percInBox(boxClassCondition.extractClass == classAssigned);
            end
            
            
            tmp_smIndxClass(iTrk,boxClassCondition.extractClass == classAssigned) = true;
            
            %tmp_maskSizeClass(iTrk, max(tmp_class)) = tmp_maskSize(find(max(tmp_class)));
            
        end % for iTrk = 1:tmp_numTrk
        
        % STEP 2: Calculate mask density of each region:
        % Calc. number of tracklets in each region ------------------------------------------------
        cntTracklet = sum(tmp_smIndxClass,1);
        
        % Calc. mask size of that region (considering only the original mask) ----------------------
        
        % Consider only the original mask (boxes that were not originally part
        % of the original mask and did not contribute to any SM detection
        % should not contribute to any mask size):
        %         ogMaskStatus = double(cell2mat(smTrackletIdPerBox{i}(iMaskFr).maskStatus) > 0); % 0 = not part of mask, 1 = part of mask regardless of edge status
        %ogMaskSize = smTrackletIdPerBox{i}(iMaskFr).maskSize;
        totMaskSize = nan(1, numClass);
        for iClass = 1:numClass
            tmp_boxEdgeClass = boxClass{i}(iMaskFr).edgeClass == boxClassCondition.extractClass(iClass);
            %totMaskSize(1,iClass) = sum(boxClass{i}(iMaskFr).maskSize(tmp_boxEdgeClass));
            totMaskSize(1,iClass) = sum(cell2mat(smTrackletIdPerBox{i}(iMaskFr).maskSize(tmp_boxEdgeClass))); % Use smTrackletIdPerBox to get only OG mask info. Compare this with boxClass.maskSize to confirm the differences.
        end % for iClass = 1:numClass
        
        % Calc. mask density density ----------------------------------------------------------------
        maskDensityByRegion = cntTracklet./totMaskSize;
        
        
        % STEP 3: Assign the SM to its correct class output --------------------------------------
        for iClass = 1:numClass
            
            smPropPerTimeIntByBoxClass{i,iClass}(iFr).smFrames = smPropPerTimeInterval{i}(iFr).smFrames;
            smPropPerTimeIntByBoxClass{i,iClass}(iFr).smOgMaskDensity = smPropPerTimeInterval{i}(iFr).maskDensity;
            smPropPerTimeIntByBoxClass{i,iClass}(iFr).smMaskClassDensity = maskDensityByRegion(1,iClass);
            smPropPerTimeIntByBoxClass{i,iClass}(iFr).percInClassAssigned = tmp_percInClassAssigned(tmp_smIndxClass(:,iClass),iClass);
            
            for iFld = 1:length(fieldList)
                if strcmp(fieldList{iFld}, 'mergeInfoSpace') || strcmp(fieldList{iFld}, 'splitInfoSpace')
                    continue % TN 20230428: skipping merge & split information for now because their info is not associated w/ SM tracklets
                else
                    smPropPerTimeIntByBoxClass{i,iClass}(iFr).(fieldList{iFld}) = tmp_sm.(fieldList{iFld})(tmp_smIndxClass(:,iClass),:);
                end
            end
            
        end % for iClass = 1:numClass
        
    end % for iFr = 1:numFr
    
end % for i = 1:numMov

end