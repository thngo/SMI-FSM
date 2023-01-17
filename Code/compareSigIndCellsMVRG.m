function [pCompareTtestR,pCompareTtestL,outlierIndx1,outlierIndx2] = compareSigIndCellsMVRG(testResults1,testResults2,pSigTtest1,pSigTtest2,spekLabelBase, smLabelBase)
%COMPARESIGINDCELLSMVRG compares multivariate regression results for two
%groups of cells, each population
%
%SYNOPSIS:  [pCompareTtestR,pCompareTtestL,outlierIndx1,outlierIndx2] = compareIndCellsMVRG(testResults1,testResults2, spekLabelBase, smLabelBase)
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
%OUTPUT    pCompareTtest: 3D array of coefficient comparison p-values
%                         between the two groups.
%                         Dimensions are (number of speckle properties) by
%                         (number of SM properties) by (number of trajectory
%                         types). Trajectory types are: immobile, confined,
%                         free, directed and ALL.
%                         Remark: this comparision is done via tailed
%                         alternative hypothesis testing
%                         (1) pCompareTtestR: "right tail" - Test against
%                         the alternative hypothesis that the population
%                         mean of x is greater than the population mean of y.
%                         (2) pCompareTtestL: "left tail" - Test against
%                         the alternative hypothesis that the population
%                         mean of x is less than the population mean of y.
%
%           outlierIndx1,outlierIndx2: outlier index as row(cell no.),
%                         column(speckle property), according to SM
%                         property(row) and trajectory(col)  
%                    
%Khuloud Jaqaman, 02/2020
%modified Tra Ngo, 09/2020
%modified Aparajita 02/2021, fit in outlier, notboxplot, and sigstar
%% Input check input
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
if nargin < 6 %input error
    error('Expected number of input is 6.');
end

%% initiate variable
%get number of cells
numCell1 = length(testResults1);
numCell2 = length(testResults2);

spVect = [spekLabelBase;spekLabelBase];% to make the same speckle props appear one after other
spekLabel = spVect(:);
smLabel = smLabelBase;
trajLabel = {'Immobile','Confined','Free','Directed','All'};

numSmLabel = length(smLabel);

%% Constant
% following x positions (xPos)jittered in a way that i can easily assign outlier indices to population1 and 2
% accomodates upto 5 speckle properties, need to expand list if more speckle properties are added
groupsSize = size(spekLabelBase,2);% the no.of x positions or speckle props (also is length(spekLabelBase))
groups = cell(1,groupsSize);% creates empty cell array corresponding to the no. of groups
% groups is a cell array containing X positions for the sigstar function
for gIndx = 1:groupsSize
    groups{1,gIndx} = [gIndx*2, gIndx*2+0.4]; % trick sigstar into placing significance over the appropriate plot.
end % to compare between population 1 and 2
xPos = cell2mat(groups)'; % to give xposition according to groups for not boxplot
groups2 = cell2mat (groups);

for gIndx2 = 1:length(groups2)
    groupsInd{1,gIndx2} = [groups2(gIndx2)-0.15, groups2(gIndx2)+0.15]; % trick sigstar into placing significance over the appropriate plot.
end


%% Output
pCompareTtestR = NaN(groupsSize,numSmLabel,length(trajLabel));
pCompareTtestL = NaN(groupsSize,numSmLabel,length(trajLabel));

%% Looping through 5 groups of trajectories in order: immobile, confined, free, directed, all together
for iTraj = 1:5
   
    %get MVRG coefficients
    if iTraj < 5 %section for different diffusion types
       
        coefMat3D1 = [];
        for iCell = 1 : numCell1
            if ~isfield(testResults1{iCell}{2}.mvrgs{iTraj},'nanFlag') && isstruct(testResults1{iCell}{2}.mvrgs{iTraj})
                coefMat3D1 = cat(3,coefMat3D1,testResults1{iCell}{2}.mvrgs{iTraj}.coef);
            else % TO-DO: missingData matrices could be made as constant at the start of compareIndCellsMVRG.m
                missingData1 = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                coefMat3D1 = cat(3,coefMat3D1,missingData1);
            end
        end
        coefMat3D2 = [];
        for iCell = 1 : numCell2
            if isstruct(testResults2{iCell}{2}.mvrgs{iTraj}) && ~isfield(testResults2{iCell}{2}.mvrgs{iTraj},'nanFlag')
                coefMat3D2 = cat(3,coefMat3D2,testResults2{iCell}{2}.mvrgs{iTraj}.coef);
            else
                missingData2 = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                coefMat3D2 = cat(3,coefMat3D2,missingData2);
            end
        end
       
    else %section for all trajectories together (regardless of diffusion type)
       
        coefMat3D1 = [];
        for iCell = 1 : numCell1
            if isstruct(testResults1{iCell}{1}.mvrgs) && ~isfield(testResults1{iCell}{1}.mvrgs,'nanFlag')
                coefMat3D1 = cat(3,coefMat3D1,testResults1{iCell}{1}.mvrgs.coef);
            else
                missingData1 = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                coefMat3D1 = cat(3,coefMat3D1,missingData1);
            end
        end
       
        coefMat3D2 = [];
        for iCell = 1 : numCell2
            if isstruct(testResults2{iCell}{1}.mvrgs) && ~isfield(testResults2{iCell}{1}.mvrgs,'nanFlag')
                coefMat3D2 = cat(3,coefMat3D2,testResults2{iCell}{1}.mvrgs.coef);
            else
                missingData2 = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                coefMat3D2 = cat(3,coefMat3D2,missingData2);
            end
        end
       
    end
   
    %compare coefficients for each single-molecule property
    for iSMProp = 1 : size(coefMat3D1,2)
       
        %collect relevant coefficients
        tmp1 = squeeze(coefMat3D1(1:end,iSMProp,:))';
        tmp2 = squeeze(coefMat3D2(1:end,iSMProp,:))';
       
        for icol = 1:groupsSize % take relevant column out and run outlier analysis-gesd
            tmp1Flag(:,icol) = isoutlier (tmp1(:,icol),'gesd');
            tmp2Flag(:,icol) = isoutlier (tmp2(:,icol),'gesd');
        end
       
        [outlierindxR1,outlierindxC1]  = find(tmp1Flag == 1);
        [outlierindxR2,outlierindxC2]  = find(tmp2Flag == 1);
        % saving outlier location temporarily as subscript
       
        outlierIndxtmp1 = find(tmp1Flag == 1);% getting linear indx of outlier, to make outliers NaN
        outlierIndxtmp2 = find(tmp2Flag == 1);
       
        outlierValue1 = tmp1(outlierIndxtmp1);% collects outliervalues for 1st population -tmp1
        outlierValue2 = tmp2(outlierIndxtmp2);% collects outliervalues for 2nd population -tmp2
       
        outlierPlotIndx1 = outlierindxC1*2;% to match with x pos ,so that scatter can plot outliers in corrcet position
        outlierPlotIndx2 = outlierindxC2*2+0.4;% jittered according to xPos
       
       
        outlierIndx1{iSMProp,iTraj} = horzcat(outlierindxR1,outlierindxC1,outlierValue1); % saving outlier index as row(cell no.),column(speckle property), according to SM property(row) and trajectory(col)
        outlierIndx2{iSMProp,iTraj} = horzcat(outlierindxR2,outlierindxC2,outlierValue2); % saving outlier index as row(cell no.),column(speckle property), according to SM property(row) and trajectory(col)
       
        % set figure limits to accomodate for future outlier and allow for more space
        yMin = min(min(tmp1(:)),min(tmp2(:))) - 0.3;
        yMax = max(max(tmp1(:)),max(tmp2(:))) + 0.3;
       
        tmp1(outlierIndxtmp1) = NaN;% will make the outliers in tmp NaN
        tmp2(outlierIndxtmp2) = NaN;% will make the outliers in tmp NaN
       
        % ttest2 for sig
        [~,pValueR] = ttest2(tmp1,tmp2,'Vartype','unequal','Tail','right');
        [pValueL]   = 1-pValueR;% pValue Left
        [pValueMin] = min(pValueR,pValueL);
        % reading from ttest -pSig
        %[pValueInd] = horzcat(pSigTtest1(:,iSMProp,iTraj)',pSigTtest2(:,iSMProp,iTraj)');
       for igrpInd  = 1:groupsSize
            pValueInd{1,igrpInd} = [pSigTtest1(igrpInd,iSMProp,iTraj),pSigTtest2(igrpInd,iSMProp,iTraj)];
       end
        pValueInd = cell2mat(pValueInd);
       
        %  automate padcat to place coefficients for the same
        %  speckle properties from the 2 different populations to
        %  be compared next to each other to be input into
        %  notboxplot whoch takes in a single matrix
       
        if groupsSize == 1 % TN added 20220207 to accomodate 1 speckle propety input
            tmp = padcat(tmp1',tmp2');
        else
            for iMix  = 1:groupsSize
                tmp{1,iMix} =  padcat(tmp1(:,iMix),tmp2(:,iMix));
            end
            tmp = cell2mat(tmp);
        end
        
        %              xPosfinal = xPos(1:length(propVect1)*2,1);% get the correct no. of x positions according to the no. of speckle props
        h = figure;
        notBoxPlot(tmp,xPos ,[],[],0)
        %                groupsfinal = groupsCompare(1,1:length(propVect1)); % which groups to plot
       
        % assign significance for sigstar , pValue as row and give groupings as cells
         sigstar(groupsInd,pValueInd)% significance indicator for individual ttests
        hold on
       
        sigstar(groups,pValueMin) % significance indicator
        hold on
       
        scatter( outlierPlotIndx1,outlierValue1)
        scatter( outlierPlotIndx2,outlierValue2)
       
        ylim([yMin yMax])
        title(['MVRG - ' smLabel{iSMProp} ' - ' trajLabel{iTraj}]);
        h.Children.XTickLabel = spekLabel;% needs to be changed
        saveas(h,['mvrg_' smLabel{iSMProp} '_' trajLabel{iTraj} '.fig']);
        close all;
       
        pCompareTtestR(:,iSMProp,iTraj) = pValueR';
        pCompareTtestL(:,iSMProp,iTraj) = pValueL';
       
        clear tmp
        clear pValueInd
    end
   
end
end