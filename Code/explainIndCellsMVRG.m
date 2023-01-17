function [pSigTtest, outlierIndx, randResult, mvrgCoefStruct] = explainIndCellsMVRG(testResults, spekLabelBase, smLabelBase, sigThreshInput)
%EXPLAININDCELLSMVRG displays multivariate regression results for individual cells and pick out outlier cells
%
%SYNOPSIS: [pSigTtest, outlierIndx, randResult, mvrgCoefStruct] = explainIndCellsMVRG(testResults, spekLabelBase, smLabelBase, sigThreshInput)
%
%INPUT     testResults  : Cell array of individual cell analysis results.
%                         This function cares about the multivariate
%                         regression output within this cell array. See
%                         smactinAggMov for details.
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
%                       5: speckle intensity (ilmax)
%                       6: speckle background intensity (ibkg)
%                       7: ilmaxCorrected
%                       8: ibkgCorrected
%                     propVect2: [vector of interger] of SM properties in tassm
%                     collected for normAggMat in testResults
%                       1: smNetSpeed
%                       2: smMeanAmp
%                       3: smDiffusionCoefficients
%                       4: smDensity
%                       5: smDiffusionRadius
%                       6: smAggregState (oligomerization state)
%
%OUTPUT  pSigTtest  : 3D array of coefficient significance p-values.
%                     Dimensions are (number of speckle properties) by
%                     (number of SM properties) by (number of trajectory
%                     types). Trajectory types are: immobile, confined,
%                     free, directed and ALL.
%        outlierindx:  Cell array of indices and outlier values for the
%                     "outliers" detected by the outlier detection algorithm "GESD".
%                      Columns are trajectory types and rows are single
%                      molecule properties.Each cell is a matrix with
%                      column1 : the movie no. column2: speckle property,
%                      column3: outlier value. Multiple rows indicate
%                      different outliers.
%       randResult : (1 x 2 cell) summary of test of significance for randomization
%                     randResult{1} = fractionFlagMat3D; 
%                     randResult{2} = avgPval; 
%       mvrgCoefStruct: (structure) containing the inliers of plotted MVRG
%                     coefficients, with the following fields:
%                           .Immobile
%                           .Confined
%                           .Free
%                           .Directed
%                           .All
%
%Khuloud Jaqaman, 02/2020
%Modify Tra Ngo, 06/2020 & 08/2021
%Modify Aparajita Dasgupta 02/2021
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

%% Input, constant, & variable initiation
if ~exist('sigThreshInput','var') || isempty(sigThreshInput)
    alphaTtest = 0.05; % alpha-value for ttest at 0.05
    alphaRandTest = 0.05; % alpha-value for randomization test at 0.05
else
    alphaTtest = sigThreshInput(1);
    alphaRandTest = sigThreshInput(2);
    if alphaTtest > 1 || alphaTtest < 0 || alphaRandTest > 1 || alphaRandTest < 0
        error('Alpha values must be between 0 and 1.')
    end
end
sigThreshRand = [alphaTtest alphaRandTest];

if nargin < 3 %input error
    error('Expected number of input is 3.');
end

if ~iscell(testResults{1})
    error('Input must contain more than 1 movie.');
end

numCell = length(testResults); % number of cells
numSmLabel = length(smLabelBase);

% groups
groupsSize = size(spekLabelBase,2); % number of x positions, this is also length(spekLabelBase)
groups = cell(1,groupsSize);% creates empty cell array corresponding to the no. of groups
% groups is a cell array containing X positions for the sigstar function
for gIndx = 1:groupsSize
    groups{1,gIndx} = [gIndx-0.15, gIndx+0.15]; % to trick sigstar into placing significance over the appropriate plot.
end

% check the dimensions of the input mvrg result against the requested label
% display. Assumption: the same sm/speckle property dependencies where
% tested for MVRG of all cells/movies from input. This check is to prevent
% mismatching inputs.
iniCell = 1;
try % get dimensions of coef matrix of free group in the first available cell/movie
    mvrgCoefDim = size(testResults{1}{2}.mvrgs{3}.coef);
catch
    while ~isstruct(testResults{iniCell}{2}.mvrgs{3}) && iniCell <= numCell
        mvrgCoefDim = size(testResults{1}{2}.mvrgs{3}.coef);
        iniCell = iniCell+1;
    end
end
if ~exist('mvrgCoefDim','var')
    error('For all cells, free motion type, no mvrg was run.');
elseif mvrgCoefDim(1) ~= groupsSize || mvrgCoefDim(2) ~= numSmLabel
    error('Dimensions of MVRG coef results is different from requested labels.');
end

%% Output initialization
trajLabel = {'Immobile','Confined','Free','Directed','All'};
pSigTtest = NaN(groupsSize,numSmLabel,length(trajLabel));
outlierIndx  = [];

% Output randResult:
randFlag = isfield(testResults{iniCell}{1},'mvrgRand'); % assumption: if 1 cell had mvrgRand, all cells would have mvrgRand.
randResult = cell(1,2); %  randResult{1} = fractionFlagMat3D; randResult{2} = avgPval;
if randFlag % 1 passes randomization; 0 does not pass.
    [~,randResult{1}, randResult{2}] = explainIndCellsMvrgWithRand(testResults, spekLabelBase, smLabelBase, sigThreshRand);
end


%% Looping through 5 groups of trajectories in order: immobile, confined, free, directed, all together
for iTraj = 1:5
    
    %get MVRG coefficients
    if iTraj < 5 %section for different diffusion types
        
        coefMat3D = [];
        for iCell = 1 : numCell
            if ~isfield(testResults{iCell}{2}.mvrgs{iTraj},'nanFlag') && isstruct(testResults{iCell}{2}.mvrgs{iTraj})
                coefMat3D = cat(3,coefMat3D,testResults{iCell}{2}.mvrgs{iTraj}.coef);
            else
                missingData = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                coefMat3D = cat(3,coefMat3D,missingData);
            end
        end
        
    else % for all trajectories together (regardless of diffusion type)
        
        coefMat3D = [];
        for iCell = 1 : numCell
            if isstruct(testResults{iCell}{1}.mvrgs) && ~isfield(testResults{iCell}{1}.mvrgs,'nanFlag')
                coefMat3D = cat(3,coefMat3D,testResults{iCell}{1}.mvrgs.coef);
            else
                missingData = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                coefMat3D = cat(3,coefMat3D,missingData);
            end
        end
        
    end
    
    %plot coefficients for each single-molecule property
    for iSMProp = 1 : size(coefMat3D,2)
        
        %collect relevant coefficients
        tmp = squeeze(coefMat3D(1:end,iSMProp,:))';
        % get yLim values for future figure
        yMin = min(tmp(:)) - 0.3; % to accomodate for future outlier and allow for more space
        yMax = max(tmp(:)) + 0.3;
        
        if isnan(yMin) || isnan(yMax)
            disp(['Motion type ' num2str(iTraj) ' SM prop ' num2str(iSMProp) ' is NaN only.'])
            continue
        end
        
        % run outlier analysis
        tmpFlag = false(size(tmp,1),groupsSize); % initialize tmpFlag
        for icol = 1:groupsSize
            
            % remove NaN before outlier detection, because sometimes isoutlier can be buggy when 10% of array > number of non-NaN data
            tmpCurrCol = tmp(:,icol);
            indxNoNan = find(~isnan(tmpCurrCol));
            
            tmpFlag(indxNoNan,icol) = isoutlier (tmpCurrCol(indxNoNan),'gesd'); % detecting outliers
            %tmpFlag(:,icol) = isoutlier (tmp(:,icol),'gesd');% detecting outliers
            
        end % ( for icol = 1:groupsSize )
        
        % for loop to take relevant column out and run outlier
        % analysis-gesd, get inliers outliers, 1 and 0s - tmpflag
        %test significance
        [outlierindxR,outlierindxC]  = find(tmpFlag == 1);
        % saving outlier location temporarily as subscript
        
        outlierIndxtmp = find(tmpFlag == 1);% getting linear indx of outlier,to makeoutliers NaN
        outlierValue = tmp(outlierIndxtmp); % collects outliervalues
        tmp(outlierIndxtmp) = NaN;% will make the outliers in tmp NaN
        
        [~,pValue] = ttest(tmp, 0, 'Alpha', alphaTtest);% pvalue for only the inliers
      
        
        %% Save outputs
        outlierIndx{iSMProp,iTraj} = horzcat(outlierindxR,outlierindxC,outlierValue); % saving outlier index as row(cell no.),column(speckle property), according to SM property(row) and trajectory(col)
        pSigTtest(:,iSMProp,iTraj) = pValue';
        mvrgCoefStruct.(trajLabel{iTraj}).(strrep(smLabelBase{iSMProp},' ','')) = tmp;
        
        %% Plot
        h = figure;
        
        notBoxPlot(tmp,[],[],[],0)% plots only inliers as rest are NaN
        
        
        %         groupsfinal = groups(1,1:length(propVect1));% which groups to plot
        %                 %  pValue as row and give groupings
        sigstar(groups,pValue) % significance indicator
        hold on
        
        scatter(outlierindxC,outlierValue) % plot outliers in a diff color, using just plot get coords and then plot)
        
        ylim([yMin yMax]) % to give space for outliers
        
        
        % Plot information of randomization:
        
        if randFlag % 1 passes randomization; 0 does not pass.
            
            randSig = randResult{1}(:,iSMProp,iTraj); % fraction of randomization (0 to 1)
                                                      % contains significant coefficients by ttest
                                                      % with alpha-value = alphaTtest.
            pValSig = randResult{2}(:,iSMProp,iTraj);
            
            tmpText = 'F = '; tmpText2 = repmat(tmpText,groupsSize,1); randSigStr = num2str(randSig,'%.3f'); text1 = horzcat(tmpText2, randSigStr);
            text(1:groupsSize, h.Children.YLim(1)*ones(1,groupsSize), text1,'HorizontalAlignment','center');
            
            space = diff(h.Children.YAxis.TickValues)/2;
            tmpText = 'p = '; tmpText2 = repmat(tmpText,groupsSize,1); pValSigStr = num2str(pValSig,'%.3f'); text2 = horzcat(tmpText2, pValSigStr);
            text(1:groupsSize, h.Children.YLim(1)*ones(1,groupsSize)-space(1), text2,'HorizontalAlignment','center');
            
            pValRealStr = num2str(pValue','%.3f'); text3 = horzcat(tmpText2, pValRealStr);
            text(1:groupsSize, h.Children.YLim(1)*ones(1,groupsSize)-space(1)*2, text3,'Color','red','HorizontalAlignment','center');
            
            yMin = yMin - space(1)*3;
            
            ylim([yMin yMax]) % to give space for info of randomized data
            
        end % ( if randFlag )
        
        title(['MVRG - ' smLabelBase{iSMProp} ' - ' trajLabel{iTraj}]);
        h.Children.XTickLabel = spekLabelBase;
        saveas(h,['mvrg_' smLabelBase{iSMProp} '_' trajLabel{iTraj} '.fig']);
        close all;
        
    end
    
    mvrgCoefStruct.speckleLabel = spekLabelBase;
    
end