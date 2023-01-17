function [thresholdsTable, inputParam, varOfResid, evalInfoBIC, biphasicThreshArr, stat, nonVshapeList, randPlotData] = ...
    scan4smSpeckleThreshold(tassm, matchindx, smactinFlag, scanInfo, motionType, figureInfo, scanParamArray)
%SCAN4SMSPECKLETHRESHOLD checks whether 1-line or 2-line fit describes the input data better
%
%SYNOPSIS:  [thresholdsTable, inputParam, varOfResid, evalInfoBIC, biphasicThreshArr, stat, nonVshapeList, randPlotData] = ...
%    scan4smSpeckleThreshold(tassm, matchindx, smactinFlag, scanInfo, motionType, figureInfo, scanParamArray)
%
%INPUT:     tassm      : (n x 1 cell array) data matrix of SM and actin properties.
%                        See smactinMat.m for details.
%
%           matchindx  : indices of speckles,ks, or sms matching to an object;
%                        empty if many neighbors matching is performed.
%                        See smactinMat.m and smactinAggMov.m for details.
%
%           smactinFlag: (structure) flags indicating methods of SMI-FSM
%                        combination. IMPORTANT: Field .randFlag should be
%                        set to 0 so that no unnecessary randomization is
%                        done within MVRG analysis. 
%                        See smactinPropPerTimeInt.m for details.
%
%           scanInfo   : (1x2 vector, double) values of thresholds to be tested
%                      scanInfo(1) : column index of dividing speckle
%                                    properties) (i.e. 3 = speckle density)
%                      scanInfo(2) : step size of scanning values (i.e. 0.0005)
%
%           motionType : (integer, optional) diffusion type to analyse. 
%                       = 0 : immobile
%                       = 1 : confined brownian
%                       = 2 : pure brownian (free diffusion)
%                             If input [], default to 2 (free).
%                       = 3 : directed motion  
%
%           figureInfo : (1 x 4 cell array) Flag indicating whether figures
%                        will be output or not.
%                       Cell {1}: (logic) 1 (Default) = output figures & save
%                                them in current directory under file
%                                names: '2lineOver1_.fig'; 'varmin_1_.fig';
%                                'varmin_2_.fig'.
%                                0 = no output figures;
%                       Cell {2}: (char) string to be included in figure's name
%                                for saving .fig into current directory.
%                       Cell {3}: (logic) 1 = output figures of 1 randomized data
%                                example & save them in current directory
%                                under file names: 'scatterRandData_.fig';
%                                0 (Default) = not output figures;
%                       Cell {4}: (logic) 1 = output figures of real data &
%                                save them in current directory under file 
%                                names: 'scatterRealData_.fig';
%                                0 (Default) = not output figures;
%
%           scanParamArray: (1x4 cell array) (FKA propStruct) Indication of 
%                           which property to perform analysis on:
%                       {1}: property column index of speckle, i.e = [3]
%                       {2}: property column index of SM, i.e = [3]
%                       {3}: max number of property: (1x2 vector) number of speckle
%                            property and number of SM property in tassm.
%                            i.e = [9 6]
%                       {4}: outlier parameters (4x1 vector), i.e  = [2 0.1 4 -1]
%                       See smactinAggMov.m for details.
%                       {5}: bootstrapping flag. Optional. Default = 0 (no
%                       	 bootstrap). If = M, then M bootstrap is
%                       	 performed on data (random sampling with
%                       	 replacement is done before performing biphasic
%                       	 analysis.) -- this option is obsolete.
%
%OUTPUT:     thresholdsTable: (table) reporting at which thresholds is fit
%                             variance at minimum, include the following column:
%                       .thresholds      : threshold at which 2-line fit has min variance
%                       .minVar          : variance of 2-line fit at listed threshold 
%                       .vShapeFlagFinal : final decision flag indicating
%                                          if 2-line fit describes the data
%                                          better than the 1-line fit (1)
%                                          or not (0).
%                       .vShapeFlagReal  : flag indicating if the 2-line
%                                          fit describes the real data
%                                          better than the 1-line fit (1)
%                                          or not (0).
%                       .vShapeFlagRand  : flag indicating if the 1-line
%                                          fit describes the randomized
%                                          data better than the 2-line fit
%                                          m < M times (1) or not (0). M is
%                                          hardcoded to 5. 
%                       .vShapeNumRand   : m from .vShapeFlagRand.
%                       .valBICvShape    : BIC value of 2-line fit.
%                       .valBICline      : BIC value of 1-line fit.
%                       .edgeFlag        : flag indicating whether the
%                                          listed threshold is at the edge
%                                          of all scanned thresholds (1) or
%                                          not (0).
%                       .coefsOfMin      : regression coefficients, each
%                                          entry is a 2x3 matrix. row 1 =
%                                          y-intercept; row 2 = slope.
%                               col1 = below threshold (2-line fit at listed threshold)
%                               col2 = above threshold (2-line fit at listed threshold)
%                               col3 = total data population (1-line fit)
%
%             inputParam     : input parameters.
%
%             varOfResid     : variance of residues from regression of both
%                              population, calculated after grouping the
%                              residues from both population.
%                              Cols: different thresholds;
%                              Row: different cells
%                              = nansum((residuals-resMean).^2)/(n-1);
%
%             evalInfoBIC    : (Cells) each contains matrix where:
%                              Rows = different movies/cell ROIs
%                              Cols = different thresholds that were scanned
%                            {1,1}: (logical matrix) indicating if 2-line
%                                   (1) describe the data better 1-line fit 
%                                   (0) based on BIC.
%                            {1,2}: (matrix) BIC value of 2-line fit
%                            {1,3}: (matrix) BIC value of 1-line fit
%
%           biphasicThreshArr: array of thresholds of movies that are
%                             deemed biphasic & not biphasic in randomized
%                             case.
%
%               stat         : (1x4 vector) statistics of biphasic behavior.
%                              1st 3 entries: [mean, median, std] of biphasicThreshArr.
%                              4th entry = percentage of cells that passed
%                              biphasic tests out of cells input in tassm.
%
%               nonVshapeList: (matrix) index(es) of cells that failed all
%                              the biphasic criteria.
%
%               randPlotData : (n x 2 array of cell) Column 1 is a nx3
%                              matrix of data plotted for the randomized
%                              scatter plot & the corresponding threshold.
%                              Column 2 is a 2x3 matrix of slopes (row2) and 
%                              intercepts (row1) of fits below the
%                              threshold (col1), above threshold (col2),
%                              and over all data (col3) (similar to
%                              .coefsOfMin in output table.)
%                              This is only output when figureInfo{3} = true.
%                              
%
%PROCESS: thresholds are automatically determined within the code so that
%there are sufficient data for 2-line fit to be performed. Fits are
%determined by visualizing residual sum of squares of regression for the
%data below & above the user-input threshold.
%
%GLOSSARY:  SMI   : single molecule imaging
%           SM    : single molecule
%           FSM   : fluorescent speckle microscopy
%           tassm : total attributes speckles and single molecules
%           BIC   : Bayesian Information Criterion
%
%EXAMPLE INPUTS:
%               smactinFlag = struct('combFlag',10,'randFlag',0, 'match',2, 'synthData',0); 
%               scanInfo = [3 0.0005];
%               motionType = 1;
%               figureInfo = {1, 'CD36',1, 1};
%               scanParamArray = {3, 3, [9 6], [2 0.1 4 -1], 10};
%
%Tra Ngo, July 2020. Modified August 2021.
%REMARK: In the future, need to load/create default matchindx and smactinFlag
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

%% Input control, Initialization & constants
if ~exist('motionType','var') || isempty(motionType)
    motionType = 2;
end

if ~exist('figureInfo','var') || isempty(figureInfo)
    figureFlag = 1;
    figIdentifier = '';
    flagRandFig = 0;
    flagRealFig = 0;
elseif length(figureInfo) == 1
    figureFlag = figureInfo{1};
    figIdentifier = '';
    flagRandFig = 0;
    flagRealFig = 0;
elseif length(figureInfo) == 2
    figureFlag = figureInfo{1};
    figIdentifier = figureInfo{2};
    flagRandFig = 0;
    flagRealFig = 0;
elseif length(figureInfo) == 3 % TN20220317
    figureFlag = figureInfo{1};
    figIdentifier = figureInfo{2};
    flagRandFig = figureInfo{3};
    flagRealFig = 0;
elseif length(figureInfo) == 4 % TN20220317
    figureFlag = figureInfo{1};
    figIdentifier = figureInfo{2};
    flagRandFig = figureInfo{3};
    flagRealFig = figureInfo{4};
else
    error('Input "figureInfo" expects either 1x2 or 1x3 cell array or empty.') % TN20210816
end

if ~exist('scanParamArray', 'var') || isempty(scanParamArray)
    error('scanParamArray of 1x5 (or 1x4) cell array must be input.')
else
    propSpek = scanParamArray{1}; % i.e. propSpek = [3]; % just speckle actin
    propSM = scanParamArray{2}; % i.e. propSM = [3]; % just SM diffusion coefficients
    maxNumSpklProp = scanParamArray{3}; % i.e. maxNumSpklProp = [9 6]; % ILmax, IBkg, ILmaxCorrected, IBkgCorrected, averageSpeed incorporated (TN20211025)
    outlierParam = scanParamArray{4}; % i.e. outlierParam = [2 0.1 4 -1]; (-1 so outlier detection is not performed here)
    
%     if length(scanParamArray) < 5
%         bootstrapNum = 0;
%     else
%         bootstrapNum = scanParamArray{5}; % i.e. bootstrapNum = 10;
%     end
    
end

if ~exist('scanInfo','var') || isempty(scanInfo) || length(scanInfo) ~= 2
    error('Input "scanValues" expects 1x2 vector: [divPropIndx, scanValSteps].')
else
    divPropIndx = scanInfo(1);
    scanValSteps = scanInfo(2);
end

nRandom = 100; % number of randomized tassm
tolRand = 5; % tolerance for number of randomization that can be vShape (if randomized data
% is vShape more than 5 times i.e. 6, that cell would not be considered a genuine vShape.)

scanValMat = [];
numCell = length(tassm);

if isempty(matchindx)
    matchindx = cell(numCell,1);
end

%% Output initialization
inputParam = {smactinFlag, scanInfo, motionType, figureFlag, scanParamArray};
thresholdsTable = array2table(zeros(0,0));

randPlotData = cell(numCell,2); % only filled in this output if flagRandFig == 1
testResultsIndx = cell(numCell,1);

%% Get scanValues ensuring each side of threshold must have 15% (const minDatPerc) of total data
minVal = NaN; maxVal = NaN;
boundMat = cell(numCell,1);
for iCell = 1:numCell
    if isempty(tassm{iCell})
        continue
    end
    motionVect = tassm{iCell}{1,2}(:,end); % SM motion type
    indxMotionType = motionVect == motionType; % indexing which row if of motion type of interest
    spkProp = tassm{iCell}{1,1}(indxMotionType, divPropIndx);
    % NOTE: (TN20220919) Density of 0, even though not counted towards MVRG
    % analysis, is still counted here as number of datapoint to determine
    % the threshold. i.e., if dataMat has 15 datapoint with density = 0,
    % then threshold would start at 0.
    
    minDatPerc = 0.15; % each side of threshold must have at least 15% of total data
    %boundMat{iCell,1} = round(quantile(spkProp,[minDatPerc (1-minDatPerc)]),5); % quantile calls prctile, also run a bit slower than prctile
    boundMat{iCell,1} = round(prctile(spkProp,[minDatPerc*100 (1-minDatPerc)*100]),5); % same output with quantile call above
    
    minVal = min(minVal, boundMat{iCell,1}(1));
    maxVal = max(maxVal, boundMat{iCell,1}(2));
end

scanValues = (minVal: scanValSteps :maxVal);
numThresh = length(scanValues);

vShapeFlag   = nan(numCell, numThresh);
valBICvShape = nan(numCell, numThresh);
valBICline   = nan(numCell, numThresh);
varOfResid = nan(numCell, numThresh);
varianceLine = nan(numCell, numThresh);
scanValMat   = nan(numCell, numThresh);
coefArr = cell(numCell, numThresh);

for iScan = 1:numThresh
    %% Threshold of consideration in this loop
    cutoffInput = [scanValues(iScan) divPropIndx];
    scanValMat(:,iScan) = ones(numCell,1) * scanValues(iScan); % for plotting
    
    for iCell = 1:numCell
        
        if isempty(boundMat{iCell,1}) || scanValues(iScan) < boundMat{iCell,1}(1) || scanValues(iScan) > boundMat{iCell,1}(2)
            % for current cell, if the scan value is outside the range to
            % search for threshold (where the majority of the data is), then
            % continue to the next cell.
            vShapeFlag(iCell, iScan)   = NaN;
            valBICvShape(iCell, iScan) = NaN;
            valBICline(iCell, iScan)   = NaN;
            varOfResid(iCell, iScan) = NaN;
            varianceLine(iCell, iScan) = NaN;
            coefArr{iCell, iScan} = NaN;
            
        else
            
            % Regression of individual cell with each threshold (from scanValues input) via smactinAggMov
            [testResultsIndx{iCell}, ~, testResultCutOff{iCell}, inPar] = ...
                smactinAggMov(tassm(iCell), matchindx(iCell),smactinFlag, ...
                propSpek, propSM, maxNumSpklProp, 1, cutoffInput, outlierParam);
            % mvrgAggMatFlag = 1 because we observe biphasic on raw data and we only regress on 1 property vs 1 other property
            
            [vShapeFlagTmp, valBICvShapeTmp, valBIClineTmp, varOfResidTmp, varLineTmp, coefTmp] = ...
                getBicAndVarResFrMvrg(testResultsIndx{iCell}, testResultCutOff{iCell}, motionType, propSpek, inPar);
            
            vShapeFlag(iCell, iScan)   = vShapeFlagTmp;
            valBICvShape(iCell, iScan) = valBICvShapeTmp;
            valBICline(iCell, iScan)   = valBIClineTmp;
            varOfResid(iCell, iScan) = varOfResidTmp;
            varianceLine(iCell, iScan) = varLineTmp;
            coefArr{iCell, iScan} = coefTmp; % TN20211025 added
            
        end % if scanValues(iScan) < boundMat{iCell,1}(1) || scanValues(iScan) > boundMat{iCell,1}(2)
        
    end % for iCell = 1:numCell
    
end % for iScan = 1:numThresh

evalInfoBIC = {vShapeFlag, valBICvShape, valBICline};

%% Get actual threshold in which variance is minimum
[minVar,minIndx] = min(varOfResid,[],2); %minVar = value of minimum of row ith, minIndx = col index of that minimum
threshOfMin = scanValues(minIndx); % threshOfMin = threshold that produce this minimum

%% Randomize tassm (data) of each cell, shuffle observations only within groups, nRamdom of times
for iCell = 1:numCell
    
    cutoffInput = [threshOfMin(iCell) divPropIndx];
    
    for iRand = 1:nRandom
        
        shuffledTassmCurr = shuffleTassmSpkl(tassm(iCell), false);
        
        [testResultsIndxRand, ~, testResultCutOffRand, inParRand] = ...
            smactinAggMov(shuffledTassmCurr, matchindx(iCell),smactinFlag, ...
            propSpek,propSM,maxNumSpklProp,1,cutoffInput, outlierParam);
        
        [vShapeFlagRand_tmp, valBICvShapeRand_tmp, valBIClineRand_tmp, ~, ~, coefVectRand_tmp] = ...
            getBicAndVarResFrMvrg(testResultsIndxRand, ...
            testResultCutOffRand, motionType, propSpek, inParRand);
        
        vShapeFlagRand(iCell, iRand) = vShapeFlagRand_tmp;
        
        % Plot randomized figure
        if flagRandFig && iRand == 1
            [xDat, yDat, ~] = plotFitLinesOverScatter(testResultsIndxRand, testResultCutOffRand, motionType, propSM, propSpek, maxNumSpklProp, threshOfMin(iCell));
            ylabel('propSmRand'); xlabel('propSpekRand');
            saveas(gcf,['scatterRandData_' figIdentifier '_' num2str(iCell) '.fig']); close all
            randPlotData{iCell,1} = padcat(xDat, yDat, threshOfMin(iCell));
            randPlotData{iCell,2} = coefVectRand_tmp; 
            randPlotData{iCell,3} = cell2struct({valBICvShapeRand_tmp,valBIClineRand_tmp}, {'BIC2line', 'BIC1line'}, 2);
        end
        
    end % for iRand = 1:nRandom
    
end % for iCell = 1:numCell

vShapeNumRandOfMin = sum(vShapeFlagRand,2);
vShapeFlagRandOfMin = vShapeNumRandOfMin > tolRand;

% Sort out the rest of results at threshold where variance of residuals is minimum.
vShapeFlagOfMin = logical.empty(numCell,0);
edgeFlag = logical.empty(numCell,0);
valBICvShapeOfMin = nan(numCell,1);
valBIClineOfMin = nan(numCell,1);
varLineOfMin = nan(numCell,1);
coefsOfMin = cell(numCell, 1);

for i = 1:numCell
    
    if isnan(vShapeFlag(i,minIndx(i))) % NaN cannot be assigned to logical vShapeFlagOfMin(i,1)
        continue
    end
    
    vShapeFlagOfMin(i,1) = vShapeFlag(i,minIndx(i));
    valBICvShapeOfMin(i,1) = valBICvShape(i,minIndx(i));
    valBIClineOfMin(i,1) = valBICline(i,minIndx(i));
    varLineOfMin(i,1) = varianceLine(i, minIndx(i)); % I could also take any col because regression of 1-line fit is the same for every threshold, but doing this for consistency
    
    % ask if from index 1 to mindIndx(i) (where min threshold is), if any of
    % the other threshold is not NaN. If any() = 0, no other threshold exist
    % from index 1:mindIndx(i), then the current min threshold is and edge
    % threshold. -TN 20210722
    if minIndx(i) == 1 || minIndx(i) == size(varOfResid,2) % if the index is 1 or size(varOfResid,2), then threshold is already at the edge of scanValues
        edgeFlag(i,1) = true;
    else
        startEdgeFlag = ~any(~isnan(varOfResid(i,1:minIndx(i)-1))); % 0 = not an edge threshold.
        endEdgeFlag = ~any(~isnan(varOfResid(i,minIndx(i)+1:end))); % 0 = not an edge threshold.
        edgeFlag(i,1) = startEdgeFlag || endEdgeFlag;
    end
    
    % Regression coefficients of the 2-line & 1-line fit (TN20211025)
    coefsOfMin{i,1} = coefArr{i,minIndx(i)};
    
end % for i = 1:numCell

thresholdsStruct.thresholds = threshOfMin';
thresholdsStruct.minVar = minVar;
thresholdsStruct.vShapeFlagFinal = vShapeFlagOfMin == 1 & vShapeFlagRandOfMin == 0;
thresholdsStruct.vShapeFlagReal = vShapeFlagOfMin;
thresholdsStruct.vShapeFlagRand = vShapeFlagRandOfMin;
thresholdsStruct.vShapeNumRand = vShapeNumRandOfMin;
thresholdsStruct.valBICvShape = valBICvShapeOfMin;
thresholdsStruct.valBICline = valBIClineOfMin;
thresholdsStruct.edgeFlag = edgeFlag;
thresholdsStruct.coefsOfMin = coefsOfMin;


%% Save output table
thresholdsTable = struct2table(thresholdsStruct);

%% Add mean, median, std of "legit biphasic" as an output and save it.
biphasicThreshArr = thresholdsTable.thresholds (thresholdsTable.vShapeFlagFinal);
stat = [median(biphasicThreshArr), mean(biphasicThreshArr), std(biphasicThreshArr), (length(biphasicThreshArr)/numCell)];
nonVshapeList = find(~thresholdsTable.vShapeFlagFinal); % cells that failed either criteria
%disp(['cells that failed either criteria: ' num2str(nonVshapeList')])

%% Scatter (real data) of each cell with threshold
if flagRealFig
    for iCell = 1:numCell
        % plotFitLinesOverScatter(testResultsIndx{iCell}, testResultCutOff{iCell}, motionType, propSM, propSpek, maxNumSpklProp, threshOfMin(iCell));
        plotBiphasicLines(testResultsIndx, motionType, [propSpek (propSM+maxNumSpklProp(1))], thresholdsTable, iCell, figIdentifier);
        ylabel('propSM'); xlabel('propSpek');
        saveas(gcf,['scatterRealData_' figIdentifier '_' num2str(iCell) '.fig']); close all
    end
end

%% Visualize RSS and conditions
if figureFlag
    
    threshOfMinVshape = threshOfMin(:,vShapeFlagOfMin);
    threshOfMinLine = threshOfMin(:,~vShapeFlagOfMin);
    
    %% At all threshold, what is the improvement of 2-line fit over 1-line fit?
    % y = dividing the variances of the 2-line fit by the variance of the the one line fit
    f1 = figure; hold on
    title('What is the improvement of 2-line fit over 1-line fit?')
    xlabel('Different cut-off thresholds');
    ylabel('Variance of residuals from 2-line fit/Variance of residuals from 1-line fit');
    
    yMat = varOfResid./varianceLine;
    % Plot vShape lines in black
    plot(scanValMat(vShapeFlagOfMin,:)', yMat(vShapeFlagOfMin,:)','k-');
    % Plot not vShape lines in grey
    plot(scanValMat(~vShapeFlagOfMin,:)', yMat(~vShapeFlagOfMin,:)','Color',[0.5 0.5 0.5],'LineStyle','-.');
    
    % Plot 'x' mark and vertical line
    xMarkVshape = minVar(vShapeFlagOfMin,:)./varLineOfMin(vShapeFlagOfMin,:);
    xMarkLine = minVar(~vShapeFlagOfMin,:)./varLineOfMin(~vShapeFlagOfMin,:);
    scatter(threshOfMinVshape, xMarkVshape,'x'); scatter(threshOfMinLine, xMarkLine,'o')
    
    % plot vShape vertical lines in red
    for i = 1:sum(vShapeFlagOfMin)
        plot([threshOfMinVshape(i), threshOfMinVshape(i)], [min(f1.Children.YLim), xMarkVshape(i)], '--','Color','r')
    end
    % plot not vShape vertical lines in grey
    for i = 1:sum(~vShapeFlagOfMin)
        plot([threshOfMinLine(i), threshOfMinLine(i)], [min(f1.Children.YLim), xMarkLine(i)], '--','Color',[0.5 0.5 0.5],'LineStyle','-.');
    end
    saveas(f1,['2lineOver1_' figIdentifier '.fig']); close all
    
    %% At which threshold is variance at minimum?
    f2 = figure; hold on
    title('At which threshold is variance (sum(y_i - y_{mean})^2)/(n-1) at minimum?')
    xlabel('Different cut-off thresholds');
    ylabel('Variance of residuals from both populations (above & below)');
    
    % Plot vShape lines in black
    plot(scanValMat(vShapeFlagOfMin,:)', varOfResid(vShapeFlagOfMin,:)','k-');
    % Plot not vShape lines in grey
    plot(scanValMat(~vShapeFlagOfMin,:)', varOfResid(~vShapeFlagOfMin,:)','Color',[0.5 0.5 0.5],'LineStyle','-.');
    
    xVshape = minVar(vShapeFlagOfMin,:);
    xLine = minVar(~vShapeFlagOfMin,:);
    scatter(threshOfMin, minVar,'x')
    % plot vShape vertical lines in red
    for i = 1:sum(vShapeFlagOfMin)
        plot([threshOfMinVshape(i), threshOfMinVshape(i)], [min(f2.Children.YLim), xVshape(i)], '--','Color','r')
    end
    % plot not vShape vertical lines in grey
    for i = 1:sum(~vShapeFlagOfMin)
        plot([threshOfMinLine(i), threshOfMinLine(i)], [min(f2.Children.YLim), xLine(i)], '--','Color',[0.5 0.5 0.5],'LineStyle','-.');
    end
    saveas(f2,['varmin_1_' figIdentifier '.fig']); close all
    
    %% At which threshold is RELATIVE variance at minimum?
    f3 = figure;
    title('At which threshold is relative variance at minimum?')
    hold on
    
    % Scale variance of the minimum variance of each cell (row) to 1 and
    % all the other variances accordingly. (Cells are not comparable to
    % each others, only comparable within its own different scanValues).
    
    % Option 0 (KJ's): divide elements of each row (each cell) by the
    % (min(variance)) of all
    minVarVal = min(varOfResid,[],'all');
    scaledVar = varOfResid./minVarVal;
    scaledMinVar = minVar./minVarVal;
    
    % Option 1: divide elements of each row (each cell) by the floor(min(variance)) of that row
    if min(minVar) > 1
        buffer = 0;
    else
        buffer = 1;
    end
    scaledVar = varOfResid./(floor(minVar) + buffer);
    scaledMinVar = minVar./(floor(minVar) + buffer);
    
    % Option 2: divide elements of each row (each cell) by the min(variance) of that row
    scaledVar = varOfResid./minVar;
    scaledMinVar = minVar./minVar;
    
    % Option 3: divide elements of each row (each cell) by the first variance of that row (normalize all lines to start at the same point)
    for i = 1:numCell
        indxFirstVar = find(~isnan(varOfResid(i,:)),1);
        if ~isempty(indxFirstVar)
            firstMinVar(i,1) = varOfResid(i,indxFirstVar);
        else
            firstMinVar(i,1) = NaN;
        end
    end
    scaledVar = varOfResid./firstMinVar;
    scaledMinVar = minVar./firstMinVar;
    
    xlabel('Different cut-off thresholds');
    ylabel('Variance of residuals from both populations (above & below)');
    % Plot vShape lines in black
    plot(scanValMat(vShapeFlagOfMin,:)', scaledVar(vShapeFlagOfMin,:)');%,'k-');
    % Plot not vShape lines in grey
    plot(scanValMat(~vShapeFlagOfMin,:)', scaledVar(~vShapeFlagOfMin,:)','Color',[0.5 0.5 0.5],'LineStyle','-.');
    
    scaledMinVarV = scaledMinVar(vShapeFlagOfMin,:);
    scaledMinVarL = scaledMinVar(~vShapeFlagOfMin,:);
    scatter(threshOfMin, scaledMinVar ,'x')
    
    % Plot vShape vertical lines in red
    for i = 1:sum(vShapeFlagOfMin)
        plot([threshOfMinVshape(i), threshOfMinVshape(i)], [0, scaledMinVarV(i)], '--','Color','r')
    end
    % Plot not vShape vertical lines in grey
    for i = 1:sum(~vShapeFlagOfMin)
        plot([threshOfMinLine(i), threshOfMinLine(i)], [0, scaledMinVarL(i)], '--','Color',[0.5 0.5 0.5],'LineStyle','-.');
    end
    
    saveas(f3,['varmin_2_' figIdentifier '.fig']); close all
    
end % (if figureFlag)

end

%% Sub-functions
function [residualsDat, rss, variance, variance2, variance3, variance4, variance5, variance6] = ...
    getRssFromMvrg(testResultMvrgStruct)
% GETRSSFROMMVRG: takes a structure of mvrg (output of mvregress from
% smactinAggMov.m) and calculate variance of residues
% Calculating variance of residues by
% (1) RSS/n			     (4) sum(yi – ymean)/n
% (2) RSS/(n-1)		 	 (5) sum(yi – ymean)/(n-1)
% (3) RSS/(n-2)			 (6) sum(yi – ymean)/(n-2)
% Remark: the residualsDat extracted .residuals field, including any NaNs.
if ~isstruct(testResultMvrgStruct) || isfield(testResultMvrgStruct,'nanFlag')
    residualsDat = [];
    rss = NaN;
    variance = NaN;
    variance2 = NaN;
    variance3 = NaN;
    variance4 = NaN;
    variance5 = NaN;
    variance6 = NaN;
    
else
    if size(testResultMvrgStruct.residuals,2) == 1
        residualsDat = testResultMvrgStruct.residuals(:,:);
    elseif size(testResultMvrgStruct.residuals,2) == 3 || size(testResultMvrgStruct.residuals,2) == 4
        residualsDat = testResultMvrgStruct.residuals(:,2); % hard-coded column 2 because biphasic behavior is observed in y = SM diff coef (which is col = 2 in mvrg.residuals)
    else
        warning('Function is not built to handle this number of sm properties.')
    end
    [rss, variance, variance2, variance3, variance4, variance5, variance6] = ...
        getVarianceFromResiduals(residualsDat);
end

% if ~isstruct(testResultCutOff{iCell, 1}{2, 1}.mvrgs{1, motionType+1}) || ...
%         isfield(testResultCutOff{iCell, 1}{2, 1}.mvrgs{1, motionType+1},'nanFlag')
%     rssAbove(iCell,iScan) = NaN;
%     varianceAbv(iCell,iScan) = NaN;
%     varianceAbv2(iCell,iScan) = NaN;
%     residualsAbv = [];
% else
%     residualsAbv = testResultCutOff{iCell, 1}{2, 1}.mvrgs{1, motionType+1}.residuals(:,2);
%     rssAbove(iCell,iScan) = nansum(residualsAbv.^2);
%     varianceAbv(iCell,iScan) = nansum((residualsAbv-nanmean(residualsAbv)).^2)/size(residualsAbv,1);
%     varianceAbv2(iCell,iScan) = rssAbove(iCell,iScan)/size(residualsAbv,1);
%
% end
%
%
% if ~isstruct(testResultCutOff{iCell, 1}{2, 2}.mvrgs{1, motionType+1})|| ...
%         isfield(testResultCutOff{iCell, 1}{2, 2}.mvrgs{1, motionType+1},'nanFlag')
%     rssBelow(iCell,iScan) = NaN;
%     varianceBel(iCell,iScan) = NaN;
%     varianceBel2(iCell,iScan) = NaN;
%     residualsBel = [];
% else
%     residualsBel = testResultCutOff{iCell, 1}{2, 2}.mvrgs{1, motionType+1}.residuals(:,2);
%     rssBelow(iCell,iScan) = nansum(residualsBel.^2);
%     varianceBel(iCell,iScan) = nansum((residualsBel-nanmean(residualsBel)).^2)/size(residualsBel,1);
%     varianceBel2(iCell,iScan) = rssBelow(iCell,iScan)/size(residualsBel,1);
% end
end

function [rss, variance, variance2, variance3, variance4, variance5, variance6] = getVarianceFromResiduals(residuals)
% Calculating variance of residues by
% (1) RSS/n			     (4) sum(yi – ymean)/n
% (2) RSS/(n-1)		 	 (5) sum(yi – ymean)/(n-1)
% (3) RSS/(n-2)			 (6) sum(yi – ymean)/(n-2)
if isempty(residuals)
    rss = NaN;
    variance  = NaN;
    variance2 = NaN;
    variance3 = NaN;
    variance4 = NaN;
    variance5 = NaN;
    variance6 = NaN;
else
    rss = nansum(residuals.^2);
    n = size(residuals(~isnan(residuals)),1);
    resMean = nanmean(residuals);
    % (1) RSS/n
    variance = rss/n;
    % (2) RSS/(n-1)
    variance2 = rss/(n-1);
    % (3) RSS/(n-2)
    variance3 = rss/(n-2);
    % (4) sum(yi – ymean)/n
    variance4 = nansum((residuals-resMean).^2)/n;
    % (5) sum(yi – ymean)/(n-1) (calculate resMean from data points so 1 less degree of freedom)
    variance5 = nansum((residuals-resMean).^2)/(n-1);
    %variance52 = var(residuals,'omitnan');
    % (6) sum(yi – ymean)/(n-2)
    variance6 = nansum((residuals-resMean).^2)/(n-2);
end
end

function [vShapeInvertFlag, coef1, coef2] = checkInvtVshapeCoef(testResultCutOffCell, motionType, inPar)
% check invert vShape: below threshold fit slope is positive, above
% threshold fit slope is negative.

if inPar{6} % inputParam{6} = mvrgAggMatFlag (if = 1, then aggMat was used for MVRG, not normalized aggMat)
    coef2 = testResultCutOffCell{2, 1}.mvrgs{1, motionType+1}.coef(2); % coef of above population
    coef1 = testResultCutOffCell{2, 2}.mvrgs{1, motionType+1}.coef(2); % coef of below population
else
    error('Function not written with running normAggMat in mind.')
end


if coef1 > 0 && coef2 < 0
    vShapeInvertFlag = true;
else
    vShapeInvertFlag = false;
end

end

function [vShapeFlag, valBICvShape, valBICline, varOfResid, variance1Line, coefVect] = getBicAndVarResFrMvrg(testResultsIndxCell, testResultCutOffCell, motionType, propSpek, inParam)
%GETBICANDVARRESFRMVRG:  take results of MVRG (testResultsInd &
%testResultCutOff) and returns information regarding BIC assessment of 1
%line fit vs vShape fit and variance of residuals from vShape fit.
%OUTPUT: coefVect = (1x2 vector), coefVect(1) = coef of regression on below
%population, coefVect(2) = coef of regression on above population
%% Check requirement:
% % min = 11 datapoints per MVRG
% All datapoints (1 line case): Make sure that MVRG ran
mvrgAllFlag = isstruct(testResultsIndxCell{2,1}.mvrgs{1,motionType+1}) && ...
    ~isfield(testResultsIndxCell{2,1}.mvrgs{1,motionType+1},'nanFlag');

% Above threshold (2 lines case)
complIndxAbv = testResultCutOffCell{2,1}.completeFlag{1,motionType+1};
inlrIndxAbv = testResultCutOffCell{2,1}.outlierFlag{1,motionType+1};
aggMatComplAbv = testResultCutOffCell{2,1}.aggMat{1,motionType+1}(complIndxAbv,:);
aggMatComplInlrAbv = aggMatComplAbv(inlrIndxAbv,:);
% Below threshold (2 lines case)
complIndxBel = testResultCutOffCell{2,2}.completeFlag{1,motionType+1};
inlrIndxBel = testResultCutOffCell{2,2}.outlierFlag{1,motionType+1};
aggMatComplBel = testResultCutOffCell{2,2}.aggMat{1,motionType+1}(complIndxBel,:);
aggMatComplInlrBel = aggMatComplBel(inlrIndxBel,:);

insuffDatFlag = min([size(aggMatComplInlrAbv,1),size(aggMatComplInlrBel,1)]) < 11;

% % Make sure that both Below & Above group has MVRG ran
mvrgAbvFlag = isstruct(testResultCutOffCell{2,1}.mvrgs{1,motionType+1}) && ...
    ~isfield(testResultCutOffCell{2,1}.mvrgs{1,motionType+1},'nanFlag');
mvrgBelFlag = isstruct(testResultCutOffCell{2,2}.mvrgs{1,motionType+1}) && ...
    ~isfield(testResultCutOffCell{2,2}.mvrgs{1,motionType+1},'nanFlag');

if insuffDatFlag || ~mvrgAbvFlag || ~mvrgBelFlag
    flagBIC = false;
else
    flagBIC = true;
end

%% Get variance of residuals of "all" (1-line fit)
if mvrgAllFlag
    residuals1Line = getRssFromMvrg(testResultsIndxCell{2,1}.mvrgs{1, motionType+1});
    [~, ~, ~, ~, ~, variance1Line] = getVarianceFromResiduals(residuals1Line);
else
    variance1Line = NaN;
end


if flagBIC
    %% Get variance of residuals for each threshold (2-line fit)
    % Each population
    residualsAbv = getRssFromMvrg(testResultCutOffCell{2, 1}.mvrgs{1, motionType+1});
    residualsBel = getRssFromMvrg(testResultCutOffCell{2, 2}.mvrgs{1, motionType+1});
    
    residualsAll  = vertcat(residualsAbv , residualsBel ); % Both population
    
    % Get variance
    [~, ~, ~, ~, ~, varOfResid] = getVarianceFromResiduals(residualsAll);
    
    %% BIC: (simple formula based on variance (& number of available datapoints) & number of parameters)
    % Discussed: number of parameters (1 line fit => 1 parameters, 2 line fits => either 2 or 3 parameters???))
    % v-shape
    numDataPtsVshape = length(residualsAll );
    residualsVshape = residualsAll ; %residualsAll
    numParamVshape = 2*(length(propSpek) + 1); % 2 lines (unnormalized data), each has 1 slope & 1 intercept in default case
    valueVshapeBIC = numDataPtsVshape*log(sum(residualsVshape.^2)/numDataPtsVshape) + log(numDataPtsVshape)*numParamVshape;
    
    % line:
    numDataPtsLine = sum(testResultsIndxCell{2, 1}.outlierFlag{1, motionType+1});
    residualsLine = getRssFromMvrg(testResultsIndxCell{2, 1}.mvrgs{1, motionType+1});
    numParamLine = length(propSpek) + 1; % 1 line, has 1 slopes & 1 intercept in default case
    valueLineBIC = numDataPtsLine*log(sum(residualsLine.^2)/numDataPtsLine) + log(numDataPtsLine)*numParamLine; % line has slope & intercept
    
    % vshape flag
    twoLineFlag  = valueVshapeBIC < valueLineBIC; % 1 if data is better described with 2-line model (v-shape)
    vShapeInvertFlag = checkInvtVshapeCoef(testResultCutOffCell, motionType, inParam); % 1 if data transition from positive slope to negative slope when value of x (speckle prop) is increased
    vShapeFlag  = twoLineFlag & vShapeInvertFlag; % 1 if data is better described with 2-line model (v-shape), 0 if data is better described with 1-line model or if 2 lines are not vShape
    valBICvShape  = valueVshapeBIC;
    valBICline  = valueLineBIC;
    
    %% Get coefficients and slope of fit line
    coef1line = testResultsIndxCell{2,1}.mvrgs{1,motionType+1}.coef; % coef1line(1) is intercept, coef1line(2) is slope
    coefAbove = testResultCutOffCell{2, 1}.mvrgs{1, motionType+1}.coef; % coefAbove(1) is intercept, coefAbove(2) is slope
    coefBelow = testResultCutOffCell{2, 2}.mvrgs{1, motionType+1}.coef; % coefBelow(1) is intercept, coefBelow(2) is slope
    coefVect = horzcat(coefBelow, coefAbove, coef1line); % row 1 are intercepts, row 2 are slopes
    
else % if number of datapoint criteria does not pass, then cell is not biphasic.
    vShapeFlag  = false;
    valBICvShape  = NaN;
    valBICline  = NaN;
    varOfResid  = NaN;
    coefVect = NaN;
end
clear flagBIC
end

function [xDat, yDat, coefsOfMin] = plotFitLinesOverScatter(testResultsInd, testResultCutOff, motionType, propSM, propSpek, maxNumSpklProp, thres)
% thres = thresholdsTable.thresholds(1);

% Data:
mvrgCell1_Above_atMinThresh = testResultCutOff{2, 1}.mvrgs{1, motionType+1};
mvrgCell1_Below_atMinThresh = testResultCutOff{2, 2}.mvrgs{1, motionType+1};
mvrgCell1All_atMinThresh = testResultsInd{2, 1}.mvrgs{1, motionType+1};

if isfield(mvrgCell1_Above_atMinThresh,'nanFlag') || isfield(mvrgCell1_Below_atMinThresh,'nanFlag') || isfield(mvrgCell1All_atMinThresh,'nanFlag')
    xDat = NaN; yDat = NaN;
    coefsOfMin = [];
    return
else
    yDat= testResultsInd{2,1}.aggMat{1,motionType+1}(:,propSM+maxNumSpklProp(1));
    xDat = testResultsInd{2,1}.aggMat{1,motionType+1}(:,propSpek);
    
    
    figure; scatter(xDat,yDat); hold on
    
    x = min(xDat):(max(xDat)-min(xDat))/100:max(xDat);
    x2lineBelow = min(xDat):(thres-min(xDat))/100:thres;
    x2lineAbove = thres:(max(xDat)-thres)/100:max(xDat);
    
    y1line = mvrgCell1All_atMinThresh.coef(1) + mvrgCell1All_atMinThresh.coef(2)*x;
    
    plot(x,y1line,'k') % threshold
    
    yBelow = mvrgCell1_Below_atMinThresh.coef(1) + mvrgCell1_Below_atMinThresh.coef(2)*x2lineBelow;
    plot(x2lineBelow,yBelow,'r')
    
    yAbove = mvrgCell1_Above_atMinThresh.coef(1) + mvrgCell1_Above_atMinThresh.coef(2)*x2lineAbove;
    plot(x2lineAbove,yAbove,'r')
    xline(thres,'--r')
    
    % Output MVRG lines
    coefsOfMin = nan(2,3);
    coefsOfMin(1,1) = mvrgCell1_Below_atMinThresh.coef(1); % below, y-intercept
    coefsOfMin(2,1) = mvrgCell1_Below_atMinThresh.coef(2); % below, slope
    coefsOfMin(1,2) = mvrgCell1_Above_atMinThresh.coef(1); % above, y-intercept
    coefsOfMin(2,2) = mvrgCell1_Above_atMinThresh.coef(2); % above, slope
    coefsOfMin(1,3) = mvrgCell1All_atMinThresh.coef(1); % all, y-intercept
    coefsOfMin(2,3) = mvrgCell1All_atMinThresh.coef(2); % all, slope

end % isfield

end