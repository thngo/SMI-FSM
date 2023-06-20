function [specklePropPerTimeIntByBoxClass, specklePropPerTimeIntByBoxNoClass] = stratifySpeckleByBoxClass(actinPropPerTimeInterval, boxClass, boxClassCondition, speckleIdPerBox)
%STRATIFYSPECKLEBYBOXCLASS: split the actin speckle data into groups based on their box classification
%
%SYNOPSIS: [specklePropPerTimeIntByBoxClass, specklePropPerTimeIntByBoxNoClass] = stratifySpeckleByBoxClass(actinPropPerTimeInterval, boxClass, boxClassCondition, speckleIdPerBox)
%
%INPUT: actinPropPerTimeInterval: (n x 1 cell array of structures) FSM arm
%                             analysis results of SMI-FSM analysis for n
%                             movies; each structure contain speckle
%                             properties for m FSM intervals. These speckle
%                             info has already been assigned into boxes.
%                             See actinPropPerTimeInterval.m and
%                             gridActinPropPerTimeInt.m for details.
%
%       boxClass: (cell array) information for each box of the grid,
%                            including classification information. See
%                            output of classifyMaskGrid.m for details.  
%
%       boxClassCondition: structure with following fields:
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
%       speckleIdPerBox: (cell array of structure) ID of actin speckle per box,
%                        number of box is dependent on imageSize and grid
%                        size. See gridActinPropPerTimeInt.m for details.
%
%OUTPUT: specklePropPerTimeIntByBoxClass: (n x m cell array) contains various
%                         properties of actin speckle, splitted into groups 
%                         based on box classification. 
%                         n = number of movie from input actinPropPerTimeInterval; 
%                         m = number of box classes user wants to extract.
%        specklePropPerTimeIntByBoxNoClass: for debugging, should be empty (TN20230501)
%
% Tra Ngo May 2023
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

numMov = size(actinPropPerTimeInterval,1);
numMask = size(boxClass,1);
numBoxId = size(speckleIdPerBox,1);
if numMov ~= numMask || numMov ~= numBoxId
    error('The movie array inputs must have similar length.')
end

numClass = length(boxClassCondition.extractClass);
specklePropPerTimeIntByBoxClass = cell(numMov,numClass);
specklePropPerTimeIntByBoxNoClass = cell(numMov, 1);

if max(boxClassCondition.extractClass) > 2
    warning('Function uses assumption of class 1 (center) & 2 (edge) for classification.')
end

for i = 1:numMov
    numFr = size(actinPropPerTimeInterval{i},1);
    
    for iFr = 1:numFr
        
        iMaskFr = iFr + boxClassCondition.firstFrame - 1;
        
        tmp_speckle = rmfield(actinPropPerTimeInterval{i}(iFr), {'kinScorePos','kinScore','speckleMaskDensity', 'indxMask'});
        fieldList = fieldnames(tmp_speckle);
        
        tmp_numSpk = length(tmp_speckle.gridPosition); % note that speckle's .gridPosition is similar to SM's .boxAssigned
        % STEP 1: For current interval, get index of each class. (Read the
        % vector indicating which box localization of SM is in and assign
        % the types to the box. Final output is the type that the SM
        % tracklet sits in majority of the time.) -------------------------
        tmp_speckleIndxClass = false(length(tmp_speckle.gridPosition), numClass);
        
        for iSpk = 1:tmp_numSpk
            
            classAssigned = boxClass{i}(iMaskFr).edgeClass(tmp_speckle.gridPosition(iSpk, 1));
            tmp_speckleIndxClass(iSpk,boxClassCondition.extractClass == classAssigned) = true;
            
        end % for iSpk = 1:tmp_numSpk
        
        % STEP 2: Calculate mask density of each region:
        % Calc. number of tracklets in each region ------------------------------------------------
        cntSpecklePerClass = sum(tmp_speckleIndxClass,1);
        
        % Calc. mask size of that region (considering only the original mask) ----------------------
                
        totMaskSize = nan(1, numClass);
        for iClass = 1:numClass
            tmp_boxEdgeClass = boxClass{i}(iMaskFr).edgeClass == boxClassCondition.extractClass(iClass);
            %totMaskSize(1,iClass) = sum(boxClass{i}(iMaskFr).maskSize(tmp_boxEdgeClass));
            totMaskSize(1,iClass) = sum(cell2mat(speckleIdPerBox{i}(iMaskFr).maskSize(tmp_boxEdgeClass))); % Use speckleIdPerBox to get only OG mask info. Compare this with boxClass.maskSize to confirm the differences.
        end % for iClass = 1:numClass
        
        % Calc. mask density density ----------------------------------------------------------------
        maskDensityByRegion = cntSpecklePerClass./totMaskSize;
        
        
        % STEP 3: Assign the SM to its correct class output --------------------------------------
        for iClass = 1:numClass
            % Put back removed fields: 'kinScorePos','kinScore','speckleMaskDensity', 'indxMask'
            specklePropPerTimeIntByBoxClass{i,iClass}(iFr).kinScorePos = actinPropPerTimeInterval{i}(iFr).kinScorePos;
            specklePropPerTimeIntByBoxClass{i,iClass}(iFr).kinScore = actinPropPerTimeInterval{i}(iFr).kinScore;
            specklePropPerTimeIntByBoxClass{i,iClass}(iFr).indxMask = actinPropPerTimeInterval{i}(iFr).indxMask;
            specklePropPerTimeIntByBoxClass{i,iClass}(iFr).spkOgMaskDensity = actinPropPerTimeInterval{i}(iFr).speckleMaskDensity;
            specklePropPerTimeIntByBoxClass{i,iClass}(iFr).spkMaskClassDensity = maskDensityByRegion(1,iClass);
            
            for iFld = 1:length(fieldList)
                specklePropPerTimeIntByBoxClass{i,iClass}(iFr).(fieldList{iFld}) = tmp_speckle.(fieldList{iFld})(tmp_speckleIndxClass(:,iClass),:);
            end
            
        end % for iClass = 1:numClass
        
        % STEP5: Debug: speckle not assigned to the mask of either
        % classification:
        % Put back removed fields: 'kinScorePos','kinScore','speckleMaskDensity', 'indxMask'
        tmp_speckleIndxNoClass = find(tmp_speckleIndxClass(:,1) == 0 & tmp_speckleIndxClass(:,2) == 0);
        for iFld = 1:length(fieldList)
            specklePropPerTimeIntByBoxNoClass{i,1}(iFr).(fieldList{iFld}) = tmp_speckle.(fieldList{iFld})(tmp_speckleIndxNoClass,:);
        end
        
        
    end % for iFr = 1:numFr
    
end % for i = 1:numMov

end