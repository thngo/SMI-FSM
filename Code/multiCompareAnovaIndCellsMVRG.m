function [sigRes] = multiCompareAnovaIndCellsMVRG(testResultCellArr,pSigTtestArray,outlierIndxArr,spekLabelBase, smLabelBase, displSpacing)
%MULTICOMPAREINDCELLSMVRG: compares multivariate regression results for multiple groups of cells by ANOVAN, sorted by motion types and SM property
%
%SYNOPSIS:  [sigRes] = multiCompareIndCellsMVRG(testResultCellArr,pSigTtestArray,outlierIndxArr,spekLabelBase, smLabelBase)
%
%INPUT     testResults1 : Cell array of individual cell analysis results for
%                         group 1. This function cares about the multivariate
%                         regression output within this cell array. See
%                         smactinAggMov for details.
%
%          testResults2 : Same as testResults1, but for group 2.
%
%          pSigTtest1   : (matrix) of p-value to show significant
%                         difference between data and normal standard.
%                         (output of explainIndCellsMVRG.m)
%
%          pSigTtest2   : same as pSigTtest1, but for group 2.
%
%          spekLabelBase: (1 x n cell) each contains name of speckle
%                         property to be plotted on the x-axis.
%          smLabelBase  : (1 x n cell) each contains name of sm property
%                         for plot title.
%        Choose appropriate speckle and single molecule properties, which
%        were used for testResults and input in an order consistent with
%        testResults.
%               Example: if you ran for spekles: [1 3 4] and for single
%               molecules: [2 3 4 6] then the following would be your input:
%               spekLabelBase = {'Speckle Speed', 'Speckle Local Density','Speckle Lifetime'};
%               smLabelBase = { 'SM Mean Amplitude','SM Diff Coef', 'SM Local Density', 'SM Aggreg State'};
%               Remark: Test results currently can have the follwing:
%                     propVect1: [vector of interger] of speckle properties in tassm
%                     collected for normAggMat in testResults
%                       1: speckleSpeed
%                       2: speckleMvmtCohere (speckle co-movement)
%                       3: speckleDensity
%                       4: speckleLifetime
%                     propVect2: [vector of interger] of SM properties in tassm
%                     collected for normAggMat in testResults
%                       1: smNetSpeed
%                       2: smMeanAmp
%                       3: smDiffusionCoefficients
%                       4: smDensity
%                       5: smDiffusionRadius
%                       6: smAggregState (oligomerization state)
%
%           displSpacing : (Optional, 1 x 2 number vector) allow user to 
%                       input  [spaceBtwEachBoxPlot, spaceBtwEachGroup] in 
%                       final output notBoxPlot figure. Default = [0.4 5]
%
%OUTPUT: multiple notBoxPlot for each motion type, each SM property will be saved in current directory.
%               
%               sigRes: (iSMProp-by-iTraj) cell array, each cell is a
%                          structure, containing following fields:
%                           .pAnovan 
%                           .tblAnovan
%                           .statsAnovanBySpek     
%                       col = iTraj in the following order: immobile,
%                             confined, free, directed motion, all together.
%
% Tra Ngo & Aparajita Dasgupta, Apr 2021 (based on compareSigIndCellsMVRG.m)
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

%% Input check input
if nargin < 5 %input error
    error('Expected number of input is 5.');
end

lenTestResult = length(testResultCellArr);
if lenTestResult > 5
    warning('Currently plots can only display up to 5 groups. Increase spacing const if display more.')
end

if ~exist('displSpacing','var') || isempty(displSpacing)
    spaceBtwEachBoxPlot = 0.4;
    spaceBtwEachGroup = 5; % change this if display more than 5 groups.
else
    spaceBtwEachBoxPlot = displSpacing(1,1);
    spaceBtwEachGroup = displSpacing(1,2);
end

%% initiate variable
sizeOfAnovaSigDisplayBar = (lenTestResult-1) * spaceBtwEachBoxPlot; % size of big sigstar bar; if only 1 group then this does not apply

%get number of cells per testResult
for iNumCell = 1:lenTestResult
    numCell(iNumCell) = length(testResultCellArr{iNumCell});
end

spVect = [spekLabelBase;spekLabelBase];% to make the same speckle props appear one after other
spekLabel = spVect(:);
smLabel = smLabelBase;
trajLabel = {'Immobile','Confined','Free','Directed','All'};

numSmLabel = length(smLabel);

%% Constant
% following x positions (xPos)jittered in a way that i can easily assign outlier indices to population1 and 2
% accomodates upto 5 speckle properties, need to expand list if more speckle properties are added
numSpeckLabel = size(spekLabelBase,2);% the no.of x positions or speckle props (also is length(spekLabelBase))

groups = cell(1,numSpeckLabel);% creates empty cell array corresponding to the no. of groups
% groups is a cell array containing X positions for the sigstar function
for gIndx = 1:numSpeckLabel
    groups{1,gIndx} = [gIndx*spaceBtwEachGroup, gIndx*spaceBtwEachGroup + sizeOfAnovaSigDisplayBar]; % trick sigstar into placing significance over the appropriate plot.
    
    % make xPos for placing each boxPlot
    xPosCell{gIndx} = [groups{gIndx}(1):spaceBtwEachBoxPlot:groups{gIndx}(2)];
end % to compare between population 1 and 2

xPos = cell2mat(xPosCell)'; % to give xposition according to groups for not boxplot


groupIdArr = [];

for i = 1:lenTestResult % aka length(groups)
    % make groupID array for anovaN test later, testResultCellArr{1} id as 1,
    % testResultCellArr{2} id as 2, etc.
    groupIdArr = vertcat(groupIdArr, ones(length(testResultCellArr{i}),1)*i);
    
end

% Position to put pSig of individual ttest for each box plot
groups2 = xPos';

for gIndx2 = 1:length(groups2)
    groupsInd{1,gIndx2} = [groups2(gIndx2)-0.15, groups2(gIndx2)+0.15]; % trick sigstar into placing significance over the appropriate plot.
end


%% Looping through 5 groups of trajectories in order: immobile, confined, free, directed, all together

for iTraj = 1:5
    for iTR = 1:lenTestResult % iTR = iTestResult
        %get MVRG coefficients
        if iTraj < 5 %section for different diffusion types
            
            coefMat3D1{iTR} = [];
            for iCell = 1 : numCell(iTR)
                if ~isfield(testResultCellArr{iTR}{iCell}{2}.mvrgs{iTraj},'nanFlag') && isstruct(testResultCellArr{iTR}{iCell}{2}.mvrgs{iTraj})
                    coefMat3D1{iTR} = cat(3,coefMat3D1{iTR},testResultCellArr{iTR}{iCell}{2}.mvrgs{iTraj}.coef);
                else % TO-DO: missingData matrices could be made as constant at the start of compareIndCellsMVRG.m
                    missingData1 = NaN(numSpeckLabel,numSmLabel);%takes care of missing data in case mvrg is not run
                    coefMat3D1{iTR} = cat(3,coefMat3D1{iTR},missingData1);
                end
            end
            
        else %section for all trajectories together (regardless of diffusion type)
            
            coefMat3D1{iTR} = [];
            for iCell = 1 : numCell(iTR)
                if isstruct(testResultCellArr{iTR}{iCell}{1}.mvrgs) && ~isfield(testResultCellArr{iTR}{iCell}{1}.mvrgs,'nanFlag')
                    coefMat3D1{iTR} = cat(3,coefMat3D1{iTR},testResultCellArr{iTR}{iCell}{1}.mvrgs.coef);
                else
                    missingData1 = NaN(numSpeckLabel,numSmLabel);%takes care of missing data in case mvrg is not run
                    coefMat3D1{iTR} = cat(3,coefMat3D1{iTR},missingData1);
                end
            end
            
            
        end %(iTraj < 5)
    end %(iTestResult)
    
    %compare coefficients for each single-molecule property
    for iSMProp = 1 : numSmLabel % size(coefMat3D1,2)
        
        for iTR = 1:lenTestResult
            
            % collect relevant coefficients
            coefTmp{iTR} = squeeze(coefMat3D1{iTR}(1:end,iSMProp,:))';
            
            % collect relevant outliers
            outlierIndxtmp{iTR} = outlierIndxArr{iTR}{iSMProp,iTraj};
            
            % For plotting outliers but they do not affect MVRG coef mean and
            % std calculations:
            [outlierSpeckleProp] =  outlierIndxtmp{iTR}(:,2); % goal: get which speckle property this outlier is
            [outlierValue] =  outlierIndxtmp{iTR}(:,3); % goal: get the MVRG coefficient values of this outlier
            
            outlierMat{iTR} = horzcat(outlierSpeckleProp, outlierValue); % outlierMat{iTR} is specific to iTraj and iSMProp (variable will be replaced in each iteration)
            
            % X-Position of outliers on notBox plot
            outlierPlotIndx{iTR} =  spaceBtwEachGroup * outlierSpeckleProp + spaceBtwEachBoxPlot*(iTR - 1); % vector of x-pos, (i.e. outlierPlotIndx1{1} = [2 2 4 6 6 6]; outlierPlotIndx1{2} = [2.4 6.4 6.4 6.4]; etc)
            
            
            for iVal = 1:length(outlierValue)
                outlierIndxBySpeckleProp = find(coefTmp{iTR} == outlierValue(iVal));% make the outliers in coefTmp{iTR} NaN
                coefTmp{iTR}(outlierIndxBySpeckleProp) = NaN;
            end
            
        end % (for iTR  = 1:lenTestResult)
        
        %% ANOVA for testing significant differences between the input groups
        coefBySpekMat = vertcat(coefTmp{:});
        for iSpeck = 1:numSpeckLabel
            [pAnovanBySpek{iSpeck}, tblAnovanBySpek{iSpeck}, statsAnovanBySpek{iSpeck}] = ...
                anovan(coefBySpekMat(:,iSpeck), {groupIdArr}, 'model', 1, 'display','off');
        end
        pAnovaMat = cell2mat(pAnovanBySpek);
        
        % Figure limits to accomodate for future outlier and allow for more space
        yMin = min(min(coefBySpekMat)) - 0.3;
        yMax = max(max(coefBySpekMat)) + 0.3;
        
        
        % reading from ttest -pSig
        
        for igrpInd  = 1:numSpeckLabel
            for iTR = 1:lenTestResult
                pValueInd(iTR,igrpInd) = pSigTtestArray{iTR}(igrpInd,iSMProp,iTraj);
            end
        end
        
        %  Make input for notBoxPlot: place coefficients for the same
        %  speckle properties from the different populations (to
        %  be compared to each other) next to each other:
        
        coefMatNBP = [];
        
        for iSpk  = 1:numSpeckLabel
            for iTR = 1:lenTestResult
                tempVect = nan(max(numCell),1);
                tempVect2 = padcat(tempVect,coefTmp{iTR}(:,iSpk)); % 1st column is useless, 2nd column is our coef that we put in coefMatNBP
                coefMatNBP = horzcat(coefMatNBP,tempVect2(:,2));
            end
        end
        
        %% Figure
        h = figure; hold on
        notBoxPlot(coefMatNBP,xPos ,[],[],0) % coefMatNBP = input for notBoxPlot, n-by-m (n = number of MVRG coefficients, m = number of speckle properties * number of input testResult groups)
        
        % assign significance for sigstar , pValue as row and give groupings as cells
        sigstar(groupsInd,pValueInd(:))% significance indicator for individual ttests
        
        sigstar(groups,pAnovaMat) % significance indicator
        
        % Plot outliers
        hold on
        for iTR = 1:lenTestResult
            scatter(outlierPlotIndx{iTR},outlierMat{iTR}(:,2)); % outlierMat{iTR}(:,2) are the values of the outliers
        end
        
        ylim([yMin yMax])
        title(['MVRG - ' smLabel{iSMProp} ' - ' trajLabel{iTraj}]);
        h.Children.XTickLabel = spekLabel;% needs to be changed
        saveas(h,['mvrg_' smLabel{iSMProp} '_' trajLabel{iTraj} '.fig']);
        close all;
        
        clear tmp
        clear pValueInd
        
        %% Output ANOVA results:
        sigRes{iSMProp,iTraj}.pAnovan = pAnovanBySpek;
        sigRes{iSMProp,iTraj}.tblAnovan = tblAnovanBySpek;
        sigRes{iSMProp,iTraj}.statsAnovanBySpek = statsAnovanBySpek;
        % TN202208005: Previous output only output the last speckle
        % properties. This was a mistake and the intention was to output
        % all of speckle properties ANOVA testing.
        
        clear pAnovanBySpek tblAnovanBySpek statsAnovanBySpek
        clear outlierMat outlierPlotIndx
        
    end % (for iSMProp = 1 : size(coefMat3D1,2))
    
end % (for iTraj = 1:5)
end