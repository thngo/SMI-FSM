function [tassm,matchindx] = smactinMat(combinePropPerTimeInt,combineConditions,smactinFlag, diffMode)
%SMACTINMAT   Converts selected fields from combinePropPerTimeInt into data matrix format
%
%SYNOPSIS: [tassm,matchindx] = smactinMat(combinePropPerTimeInt,combineConditions,smactinFlag, diffMode)
%
%INPUT: combinePropPerTimeInt: (n x 1 cell array of structures) various properties
%                              of both single molecules and actin. 
%                              (see smactinPropPerTimeInt.m for details.)
%
%         combineConditions  : (structure) parameters for combinding SM
%                              tracks from SMI and speckles from FSM.
%                              (see smactinPropPerTimeInt.m for details.)
%
%               smactinFlag  : (structure) flags indicating methods of SMI-FSM combination 
%                              (see smactinPropPerTimeInt.m for details.)
%                              Contains the following field:
%                              .combFlag: (vector of integers from 1 to 10) flag
%                                  indicating how to combine SMI-FSM properties 
%                                  IMPORTANT: Only .combFlag = 10 is
%                                  updated.
%       
%               diffMode     : (logical) Flag indicating whether diffusion
%                              mode analysis was performed upstream of
%                              SMI-FSM analysis. Input diffMode is only
%                              implemented in cases when
%                              smactinFlag.combFlag(iMov)={7,10}
%                              Optional. Default = 0, indicating that
%                              diffusion mode analysis was NOT performed on
%                              tracks in smPropPerTimeInt.m.
%                              = 1 indicates that diffModeAnalysis was performed.
%
%OUTPUT: tassm :   (n x 1 cell array of cells) Each cell contains 3 x 1 
%                  cells, containing data matrix of SM tracklets and 
%                  speckle properties:
%                       Cell {1,1}: (k observations x 9 properties) Data 
%                       matrix of speckle properties with columns as:
%                       1: Local speckle displacement magnitude (also 'speckle speed')
%                       2: Speckle co-movement (also 'speckle movement coherence)
%                       3: Local speckle density 
%                               > expect 0 if SM is not matched to speckle
%                               > expect NaN if density value is outlier
%                       4: Speckle lifetime
%                       5: Local speckle intensity (also 'ILmax')
%                       6: Local speckle background intensity (also 'IBkg')
%                       7: speckle corrected intensity (also 'ILmaxCorrected')
%                       8: speckle corrected background intensity (also 'IBkgCorrected')
%                       9: Local speckle average displacement magnitude
%
%                       Cell {2,1}: (k observations x 7 properties) Data 
%                       matrix of SM properties with columns as:
%                       1: SM net speed
%                       2: SM intensity
%                       3: SM diffusion coefficient
%                       4: SM density
%                       5: SM spatial span (previously called "smDiffusionRadius")
%                       6: SM apparent assembly state (oligomerization state)
%                       7: SM diffusion type
%
%                       Cell {3,1}: (k observation x 2) SM tracklet mean
%                       position as reported by field .smMeanPos: (x/y-coordinates) 
%
%      matchindx: indices of speckles,ks, or sms matching to an object;
%                         empty if many neighbors matching is performed.
%
%GLOSSARY:  SMI   : single molecule imaging
%           SM    : single molecule
%           FSM   : fluorescent speckle microscopy
%           tassm : total attributes speckles and single molecules
%
%REMARKS: Properties available for concatenation are (what gets concatenated 
%                     depends on combFlag:
%
%              kinetic score                kinetic score density
%
%              speckle speed                speckle density
%              speckle movement coherence
%
%              sm mean displacement              sm mean amplitude
%              sm diffusion coefficient          sm confinement radius (applicable only for immobile/confined sms)
%              sm density (not used currently)   sm net speed   
%              sm comovement angle (not used currently)
%
% Deryl Tschoerner, May 2018. Modified by Tra Ngo & Khuloud Jaqaman.
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

%% Check input and generate outputs
if ~exist('diffMode ','var')|| isempty(diffMode)
    diffMode = 0;
end

numMovs = size(combinePropPerTimeInt,1);
tassm = cell(numMovs,1);
matchindx = tassm; % this variable is only used when smactinFlag.match == 1 

    % TN20220728: remove obsolete option smactinFlag.match == 1
    %
    % if smactinFlag.match == 1    
    %     iMatch = 1;
    % elseif smactinFlag.match == 2
    %     iMatch = 2;
    % else
    %     error('error:: invalid match parameter')
    % end

iMatch = 2;

%change value depending on the combineProperties function as input
for iMov = 1 : numMovs
    
    wind1st = combineConditions.firstFsmFrame(iMov);
    if isempty(combinePropPerTimeInt{iMov})
        continue
    end

    for iComb = 1 : length(smactinFlag.combFlag(iMov)) %1 : size(combinePropPerTimeInt,2)
        
        switch smactinFlag.combFlag(iMov)
%             case 3
%                 %concatenate nn: sm, object: speckle attributes
%                 tassm{iMov,iComb}{1} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanDisp) ... (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed) ... (:,2)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp) ...  (:,3)
%                     %vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smDensity) ... (:,4)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef) ...   (:,5)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad)...     (:,6)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];  %(:,7)
%                 
%                 tassm{iMov,iComb}{2} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed) ...     (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).matchAngle) ...       (:,2)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity)]; ... (:,3)
%                     
%             case 6
%                 %concatenate nn: sm, object: ks attributes
%                 tassm{iMov,iComb}{1} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanDisp) ... (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed) ... (:,2)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp) ...  (:,3)
%                     %vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).smDensity) ...    (:,4)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef) ...   (:,5)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad)...     (:,6)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];  %(:,7)
%                 
%                 tassm{iMov,iComb}{2} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore) ...       (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScoreDensity)]; %(:,2)
%                 
%             case 2
%                 %concatenate nn: ks, object: actin attributes
%                 tassm{iMov,iComb}{1} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore) ...       (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScoreDensity)]; %(:,2)
%                 
%                 tassm{iMov,iComb}{2} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed) ...     (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).matchAngle) ...       (:,2)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity)]; ... (:,3)
%                     
%             case 8
%                 %concatenate nn: sm, object: ks attributes
%                 tassm{iMov,iComb}{1} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore)]; ...       (:,1)
%                     ...vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).kinScoreDensity)];  %(:,2)
%                     
%                 tassm{iMov,iComb}{2} = [ ...
%                     ...vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).smMeanDisp) ... (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed) ... (:,2)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp) ...  (:,3)
%                     ...vertcat(combinePropPerTimeInt{iMov,iComb}{iNN}(wind1st : end).smDensity) ...  (:,4)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef) ...   (:,5)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad)...     (:,6)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];  %(:,7)
%                 
%             case 4
%                 %concatenate nn: actin, object: ks attributes
%                 tassm{iMov,iComb}{1} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed) ...  (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).matchAngle) ...    (:,2)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity)]; %(:,3)
%                 
%                 tassm{iMov,iComb}{2} = [ ...
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScore) ...       (:,1)
%                     vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).kinScoreDensity)]; %(:,2)
                
            case {7,10}
                %concatenate nn: actin, object: sm attributes
                tassm{iMov,iComb}{1} = [ ...
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeed), ...         (:,1)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleMvmtCohere), ...    (:,2)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleDensity), ...      %(:,3)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleLifetime), ...     %(:,4)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).ILmax), ...               %(:,5)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).IBkg), ...                %(:,6)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).ILmaxCorrected), ...      %(:,7)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).IBkgCorrected),...        %(:,8)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).speckleSpeedAveMag)];     %(:,9)
                
                if diffMode %if diffusion mode was run on smPropPerTimeInt then output tassm should have diffCoef and diffRad from diffusion mode analysis
                    diffCoefVect = vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smDiffCoefF2F);
                    diffRadVect  = vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smDiffRad);
                else %if MSS diffusion analysis was run on smPropPerTimeInt then output tassm should have diffCoef and diffRad from MSS diffusion analysis
                    diffCoefVect = vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).diffCoef);
                    diffRadVect  = vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).confRad);
                end
                
                tassm{iMov,iComb}{2} = [ ...
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smNetSpeed), ...    (:,1)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanAmp), ...     (:,2)
                    diffCoefVect, ...                                                                    (:,3)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smDensity), ...     (:,4)
                    diffRadVect, ...                                                                     (:,5)
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smAggregState) ];  %(:,6)
                
                % REMEMBER: trackClass always have to be the last column!
                tassm{iMov,iComb}{2} = [tassm{iMov,iComb}{2} ...
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).trackClass)];         % last column!
                
                tassm{iMov,iComb}{3} = [ ...
                    vertcat(combinePropPerTimeInt{iMov,iComb}{iMatch}(wind1st : end).smMeanPos) ];
                
            %{
            % not sure what's going on here?? - TN20201013
            case 10
                %concatenate nn: actin/ks, object: sm attributes
                tassm{iMov,iComb}{1} = [tassm{iMov,8}{1} tassm{iMov,7}{1}]; ... (:,1:3)
                    tassm{iMov,iComb}{2} = tassm{iMov,8}{2};                       %(:,1:5)
            %}
        end
        
        if iMatch == 1
            
            %during nn analysis, bisect nns that match or unmatch under
            %combineConditions. otherwise this should not be processed
            if iComb == 10
                
                matchindx{iMov,iComb} = matchindx{iMov,7} & matchindx{iMov,8};
            else
                %%% PUT BACK THE iMATCH %%%
                matchindx{iMov,iComb} = vertcat(combinePropPerTimeInt{iMov,iComb}(wind1st : end).dist_match) <= ...
                    combineConditions.matchRadius(iMov);
            end
        end
    end
end

end