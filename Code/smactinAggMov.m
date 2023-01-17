function [testResults,trajClass, testResults_actinCutOff, inputParam] = ...
    smactinAggMov(tassm, matchindx, smactinFlag, propVect1, propVect2, ...
    maxNumPropVect, mvrgAggMatFlag, cutOffInput, outlierParam)
%SMACTINAGGMOV        Data analysis of multiple single-molecule and actin movies
%
%SYNOPSIS: [testResults,trajClass, testResults_actinCutOff, inputParam] =
%                     smactinAggMov(tassm,matchindx,smactinFlag, propVect1, ...
%                     propVect2, maxNumSpklProp, mvrgAggMatFlag, actinCutOff, outlierParam)
%
%INPUT:        tassm: (n x 1 cell array of cells) data matrix of SM and
%                     speckle properties from n movies. 
%                     See smactinMat.m for details.
%
%          matchindx: indices of speckles,ks, or sms matching to an object;
%                     matchindx will be entirely empty if smactinFlag.match(iMov) == 2,
%                     or many neighbors matching is performed. As of
%                     20220728, smactinFlag.match(iMov) is considered
%                     removed and any code that required this input has 2
%                     hardcoded. Otherwise, matchindx contains indices of
%                     sm/actin properties if nearest neighbor analysis was
%                     performed (movie specific); 
%
%        smactinFlag: structure of flag with various fields imposing
%                     restrictions on the computation (see smactinPropPerTimeInt)
%                     .combFlag: indicates how sm and speckle are matched
%                     .synthData: indicates whether the data was
%                       synthetic or not. DISCLAIMER1: As of July 2019, we
%                       focus on the case of synthData = 0, so any sections
%                       related to synthData = 0 may be more
%                       updated/detailed; any sections related to synthData
%                       = 1 may be obsolete. 
%                     .randFlag: indicates whether MVRG on randomized data
%                       should be performed 100 times (TN20210819).
%                       randFlag = 0 means no randomization is done.
%                       randFlag = n (positive integer) means randomization is done n times.
%
%           propVect1: (vector of interger) of speckle properties in tassm collected
%                      for normAggMat. (optional) If input as empty [] then default to [1 3 4]
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
%           propVect2: (vector of interger) of SM properties in tassm collected 
%                     for normAggMat. (optional) If input as empty [] then default to [2 3 4]
%                       1: SM net speed
%                       2: SM intensity
%                       3: SM diffusion coefficient
%                       4: SM density
%                       5: SM spatial span (previously called "smDiffusionRadius")
%                       6: SM apparent assembly state (oligomerization state)
%                       7: SM diffusion type
%
%           maxNumPropVect: (1x2 vector) number of speckle property in tassm
%                         and number of SM property in tassm. 
%                         Example: maxNumPropVect = [9 6]; % 9 speckle props, 6 sm props
%
%           mvrgAggMatFlag: (logic, optional) flag indicating whether
%                           normalization should be performed before MVRG.
%                           default. 0 = MVRG on normAggMat. (normalized data)
%                                    1 = MVRG on aggMat. (original data;
%                                      cutting through intercept [oneVect X-data])
%
%           cutOffInput: (1x2 vector, optional) threshold at which data is 
%                        split into 2 subsets before any MVRG and
%                        statistical analysis.
%                           (default) empty [] -> data not split into 2
%                           1st entry: threshold dividing the data into 2.
%                           2nd entry: the column in .aggMat that the threshold 
%                           is scanned within, to split the data into 2.
%                           Example: cutOffInput = [0.01 3];
%
%           outlierParam: (4 x 1 vector, optional) containing in order 
%                         [neighborNum, iterThres, maxIter, k], where: 
%                   neighborNum = (integer) number of neighbor used for
%                                 outlier detection (default = 2).
%                     iterThres = (double) tolerance threshold for deciding 
%                                 whether another iteration of outlier detection 
%                                 should be done (default = 0.1 = 10%).
%                       maxIter = (integer) maximum number of iteration for
%                                 outlier detection. (default = 4). If
%                                 maxIter = 1, then outlier detection is 
%                                 done only once.
%                             k = (integer) Roughly, for a certain value of k, 
%                                 observations that are k*sigma away from
%                                 the mean will be considered outliers.
%                                 (default = 5, determined by default input
%                                 of functions called by subfunction 
%                                 processCompltOutliersFrMat). If k = -1,
%                                 no outlier detection is performed.
%                   Default: [2 0.1 4 -1] (empirically determined)
%           
%OUTPUT: testResults : (n x 1 cells array) each cell contains structures
%                      with the following fields of different data analysis 
%                      results for input data matrices from individual cells.
%
%         .pca       : results of principal component analysis. Contains
%                      PCs and percentage of total variance explained by each PC.
%
%         .mvrgs     : results of MVRG; contains general linear model
%                      outputs in the following fields:
%                          .coef         : general linear model coefficients
%                          .covarResponse: covariance matrix of response variables
%                          .residuals    : matrix of residuals for the fit
%                          .covarParam   : covariance matrix of MVRG coefficients
%                      Note that within the code, max iterations is
%                      currently hardcoded to be 100 and can be changed for
%                      the estimation algorithm used.
%
%         .normAggMat: contains normalized data values for different
%                      properties (defined by propVect1 and propVect2).
%                      Each row is an SM-speckle observation pair.
%
%         .aggMat    : aggregated matrix contains all properties (col);
%                      Each row is an SM-speckle observation pair. aggMat 
%                      does not contain any trajClass information.
%
%           trajClass: (n x 1 cells array) each cell contains array
%                      corresponding to each SM track's diffusion type.
%                       = 0 : immobile
%                       = 1 : confined brownian
%                       = 2 : pure brownian (free diffusion)
%                       = 3 : directed motion  
%                       = NaN: not classified.
%
%       testResults_actinCutOff: (output only available if actinCutOff is inputted)
%                                (default is empty [] if no actinCutOff is inputted)
%                                cells containing structures with fields of different
%                                testing results. 1st col = above threshold
%                                inputted; 2nd col = below threshold.
%
%       inputParam : cell of input parameters for note keeping
%
%GLOSSARY:  SMI   : single molecule imaging
%           SM    : single molecule
%           FSM   : fluorescent speckle microscopy
%           tassm : total attributes speckles and single molecules
%           MVRG  : multivariate regression
%
%Deryl Tschoerner, May 2018
%
%DISCLAIMER2: If object/nn (i.e. speckles or sm) values do not exist (i.e.
%NaN), this smactinAggMov function will run but output will be empty;
%because this function is looking at contribution of speckles and sm
%together. No meaningful information will come out if you don't have either
%of them.
%
% Tra Ngo, July 2019/May 2020
% Remark: need to set up error message if number of sm attributes is unequal in at least 1 movie
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

%% Output
testResults = cell(1,2); %[]; %{1,1} = []; testResults{1,2} = [];
trajClass = [];
testResults_actinCutOff = cell(2,2); %[];
%% Input & Constant
if smactinFlag.randFlag < 1
    randMvrgInfo = {false, 0}; % default to no randomization is done.
else 
    randMvrgInfo = {true, round(smactinFlag.randFlag)}; % randomize sm-speckle matching and perform mvrg on them n times (n = smactinFlag.randFlag)
end

if ~exist('maxNumPropVect','var')|| isempty(maxNumPropVect) || size(maxNumPropVect,2) ~=2
    error('User must input correct maximum number of speckle & SM property listed in input tassm. Expect 1x2 vector.');
    % Assumption: in all cells in tassm, number of speckle properties being considered is the same
else
    numTotProp = sum(maxNumPropVect);
    maxNumSpklProp = maxNumPropVect(1);
end

if ~exist('mvrgAggMatFlag','var')|| isempty(mvrgAggMatFlag)
    mvrgAggMatFlag = 0; % 0 = MVRG on normAggMat, 1 = MVRG on aggMat
end

if mvrgAggMatFlag ~= 0 && mvrgAggMatFlag ~= 1
    error("Invalid mvrgAggMatFlag input. Please enter only 0 or 1.");
end

if ~exist('cutOffInput','var')|| isempty(cutOffInput)
    actinCutOff = []; % actinCutOff is empty, data not split into 2
    dividingDatFlag = false;
    % since actinCutOffFlag is false, variable "dividingSpklCol" will not
    % be called anywhere.
elseif length(cutOffInput) == 1
    error('cutOffInput requires info of (1) threshold and (2) which property to divide data by.')
else
    actinCutOff = cutOffInput(1);
    dividingSpklCol = cutOffInput(2);
    dividingDatFlag = true;
end

if ~exist('propVect1','var') || isempty(propVect1)
    propVect1 = [1 3 4];
end

if ~exist('propVect2','var') || isempty(propVect2)
    propVect2 = [2 3 4];
end

% check if user input diffusion coefficient (propVect2 = 3) - added TN20210609
[diffCoefFlag, diffCoefIndx] = ismember(3,propVect2); % Assumption: speckle density is on column 3.
if  diffCoefFlag % if diffusion coefficient is input:
    % remove identical datapoints ("twins") in .normAggMat for outlier detections
    indxColRmv = length(propVect1) + diffCoefIndx; 
else
    indxColRmv = 0;
end


propVect = horzcat(propVect1 ,propVect2 + maxNumSpklProp);
%disp(strcat("Performing analysis for speckle properties #", num2str(propVect1)));
%disp(strcat("Performing analysis for single-molecule properties #", num2str(propVect2)));

numSpklProp = length(propVect1);

% input for outlier detection:
if ~exist('outlierParam','var') || isempty(outlierParam)
    outlierParam = [2 0.1 4 -1];
elseif length(outlierParam) ~= 4
    error('Unexpected parameters for outliers. Vector of size 4 expected.')
end
neighborNum = outlierParam(1); % constant for outlier detection
iterThres = outlierParam(2); % 10%
maxIter = outlierParam(3); % constant for outlier detection
k = outlierParam(4); % constant for outlier detection

%% InputParam
inputParam = {matchindx,smactinFlag, propVect1, propVect2, maxNumSpklProp, mvrgAggMatFlag, cutOffInput, outlierParam};

%% Initiate variables
numMovs = size(tassm,1);
aggMat = cell(size(tassm,2),1);
minNumData = 10; % Each MVRG has to have at least (10+1) datapoint to run
smMeanPos = [];

if length(unique(smactinFlag.combFlag)) == 1
    iComb = 1;
else
    error('Function only focuses on 1 combination, smactinFlag.combFlag == 7 (or 10), at a time. If multiple combFlag is used, consider changing code to loop through all possible combinations in movie channel.')
end

nnMat = [];
obMat = [];

%% Loop through all movies & aggregate properties into aggMat
for iMov = 1 : numMovs
    
    % counting index of all sm in a movie
    if ~isempty(tassm{iMov,iComb}) % && smactinFlag.match(iMov) == 2
        matchindx{iMov,iComb} = 1 : size(tassm{iMov,iComb}{2},1);
    elseif isempty(tassm{iMov,iComb})
        continue
    end
    
    % All current properties of speckle (nn) and sm (object)
    spklPropCurr = tassm{iMov,iComb}{1}(matchindx{iMov,iComb},:);
    smPropCurr =   tassm{iMov,iComb}{2}(matchindx{iMov,iComb},:);
    try
        smMeanPosCurr= tassm{iMov,iComb}{3}(matchindx{iMov,iComb},:);
    catch
        smMeanPosCurr = [];
    end
    % Get trajectory classes matched nn (Assumption: in smactinMat.m, trajClass is appended to the last col.)
    trajClassCurr = smPropCurr(:,end);
    trajClass = vertcat(trajClass,trajClassCurr);
    
    if smactinFlag.synthData(iMov) == 1
        warning('Starting July 2019, development of function only focus on case smactinFlag.synthData == 0.')
        aggMat{iComb} = vertcat(aggMat{iComb}, ...
            horzcat( vertcat(nnMat,spklPropCurr), vertcat(obMat,smPropCurr) )  );
    elseif smactinFlag.synthData(iMov) == 0
        % For a particular combination, loop through all movies to get matched nn/object properties
        aggMat{iComb} = vertcat( aggMat{iComb} , ...
            horzcat( vertcat(nnMat,spklPropCurr), vertcat(obMat,smPropCurr(:,1 : end - 1)) ));
    end
    % at this point, all AggMat available as aggMat{iComb}
    smMeanPos = vertcat(smMeanPos, smMeanPosCurr);
    
end % (for iMov = 1 : numMovs)
testResults{iComb, 1}.smMeanPos = smMeanPos;

%% Make density of unmatched receptor-speckle NaN, get completeLabel & outlierLabel, get normalized aggMat
[dataMatOut, completeFlag, normAggMatInlier, outlierLabelVect, normDataComplete] = ...
    processCompltOutliersFrMat(aggMat{iComb},propVect, neighborNum, iterThres, maxIter, k, indxColRmv,numTotProp);
testResults{iComb,1}.completeFlag = completeFlag;
testResults{iComb,1}.outlierFlag = outlierLabelVect;
testResults{iComb,1}.normAggMatTemp = normDataComplete;
testResults{iComb,1}.normAggMat = normAggMatInlier;
testResults{iComb,1}.aggMat = dataMatOut;

%% Divide aggMat, smMeanPos, completeFlag into different motion types ("trajClass");
%% Normalized each matched properties within the subset of each motion type
%% Detect outliers within the subset of each motion type

%% Output results in testResult
for iMotType = 0 : 3
    indxType = trajClass == iMotType;
    aggMatCurr = aggMat{iComb}(indxType,:);
    try
        meanPosCurr = smMeanPos(indxType,:);
    catch
        meanPosCurr = [];
    end
    
    testResults{iComb,2}.smMeanPos{iMotType + 1} = meanPosCurr;
    
    [dataMatOutTemp, completeFlagTemp, normAggMatInlierTemp, outlierLabelVectTemp, normDataCompleteTemp] = ...
        processCompltOutliersFrMat(aggMatCurr,propVect, neighborNum, iterThres, maxIter, k,indxColRmv,numTotProp);
    testResults{iComb,2}.aggMat{iMotType + 1} = dataMatOutTemp;
    testResults{iComb,2}.completeFlag{iMotType + 1} = completeFlagTemp;
    testResults{iComb,2}.normAggMat{iMotType + 1} = normAggMatInlierTemp;
    testResults{iComb,2}.outlierFlag{iMotType + 1} = outlierLabelVectTemp;
    testResults{iComb,2}.normAggMatTemp{iMotType + 1} = normDataCompleteTemp;
    
%     % Counting for quality control during code development:
%     disp(['Motion Type #' num2str(iMotType) '**** Num observation || Num complete observation || Num inliers']);
%     try
%         perIn = size(normAggMatInlierTemp,1)*100/size(normDataCompleteTemp,1);
%         disp([num2str(size(aggMatCurr,1)) ' || ' num2str(size(normDataCompleteTemp,1)) ' || ' num2str(size(normAggMatInlierTemp,1)) '(' num2str(perIn) '%)']);
%     catch
%     end
    
end % (for iMotType = 0 : 3)

%% Split preprocessed (complete & inlier) dataset into 2 populations by speckle density (aggMat col 3) if needed
if dividingDatFlag
    % take all data, get inlier & complete data, split that into 2
    [abvAggMat, abvPos, belAggMat, belPos] = ...
        getCompltInlierIndex(testResults{iComb,1}.aggMat, ...
        testResults{iComb,1}.completeFlag, testResults{iComb,1}.outlierFlag, ...
        testResults{iComb,1}.smMeanPos, actinCutOff, dividingSpklCol,numTotProp);
    
    % note for if which column is above cutoff ever become ambiguous
    %testResults_actinCutOff{1,1}.note = 'Above';
    
    testResults_actinCutOff{1,1}.smMeanPos = abvPos;
    testResults_actinCutOff{2,1}.smMeanPos = belPos;
    testResults_actinCutOff{1,1}.aggMat = abvAggMat;
    testResults_actinCutOff{2,1}.aggMat = belAggMat;

    %%%%%%%%%%%%%%%%% above cut-off
    [~, completeFlagTemp, normAggMatInlierTemp, outlierLabelVectTemp, normDataCompleteTemp] = ...
        processCompltOutliersFrMat(testResults_actinCutOff{1,1}.aggMat,propVect, neighborNum, ...
        iterThres, maxIter, -1, indxColRmv, numTotProp); % k = -1 so that no more outlier detection is run
    testResults_actinCutOff{1,1}.completeFlag = completeFlagTemp;
    testResults_actinCutOff{1,1}.normAggMat = normAggMatInlierTemp;
    testResults_actinCutOff{1,1}.outlierFlag = outlierLabelVectTemp;
    testResults_actinCutOff{1,1}.normAggMatTemp = normDataCompleteTemp;
    
    %%%%%%%%%%%%%%%%% below cut-off
    [~, completeFlagTemp, normAggMatInlierTemp, outlierLabelVectTemp, normDataCompleteTemp] = ...
        processCompltOutliersFrMat(testResults_actinCutOff{2,1}.aggMat,propVect, neighborNum, ...
        iterThres, maxIter, -1, indxColRmv, numTotProp); % k = -1 so that no more outlier detection is run
    testResults_actinCutOff{2,1}.completeFlag = completeFlagTemp;
    testResults_actinCutOff{2,1}.normAggMat = normAggMatInlierTemp;
    testResults_actinCutOff{2,1}.outlierFlag = outlierLabelVectTemp;
    testResults_actinCutOff{2,1}.normAggMatTemp = normDataCompleteTemp;
    
    % Splitting data of each motion types: (i = motion type)
    for i = 1:length(testResults{iComb,2}.aggMat) 
        
        [abvAggMat, abvPos, belAggMat, belPos] = ...
            getCompltInlierIndex(testResults{iComb,2}.aggMat{1,i}, ...
            testResults{iComb,2}.completeFlag{1,i}, testResults{iComb,2}.outlierFlag{1,i}, ...
            testResults{iComb,2}.smMeanPos{1,i}, actinCutOff, dividingSpklCol, numTotProp);
        % Note for if which column is above cutoff ever become ambiguous:
        % testResults_actinCutOff{1,2}.note{1,i} = 'above';
        testResults_actinCutOff{1,2}.smMeanPos{1,i} = abvPos;
        testResults_actinCutOff{2,2}.smMeanPos{1,i} = belPos;
        testResults_actinCutOff{1,2}.aggMat{1,i} = abvAggMat;
        testResults_actinCutOff{2,2}.aggMat{1,i} = belAggMat;
        
        %%%%%%%%%%%%%%%%% above cut-off
        [~, completeFlagTemp, normAggMatInlierTemp, outlierLabelVectTemp, normDataCompleteTemp] = ...
            processCompltOutliersFrMat(testResults_actinCutOff{1,2}.aggMat{1,i} ,propVect, neighborNum, ...
            iterThres, maxIter, -1, indxColRmv, numTotProp); % k = -1 so that no more outlier detection is run
        testResults_actinCutOff{1,2}.completeFlag{1,i} = completeFlagTemp;
        testResults_actinCutOff{1,2}.normAggMat{1,i} = normAggMatInlierTemp;
        testResults_actinCutOff{1,2}.outlierFlag{1,i} = outlierLabelVectTemp;
        testResults_actinCutOff{1,2}.normAggMatTemp{1,i} = normDataCompleteTemp;
        
        %%%%%%%%%%%%%%%%% below cut-off
        [~, completeFlagTemp, normAggMatInlierTemp, outlierLabelVectTemp, normDataCompleteTemp] = ...
            processCompltOutliersFrMat(testResults_actinCutOff{2,2}.aggMat{1,i} ,propVect, neighborNum, ...
            iterThres, maxIter, -1, indxColRmv, numTotProp); % k = -1 so that no more outlier detection is run
        testResults_actinCutOff{2,2}.completeFlag{1,i} = completeFlagTemp;
        testResults_actinCutOff{2,2}.normAggMat{1,i} = normAggMatInlierTemp;
        testResults_actinCutOff{2,2}.outlierFlag{1,i} = outlierLabelVectTemp;
        testResults_actinCutOff{2,2}.normAggMatTemp{1,i} = normDataCompleteTemp;
    end
end

%% =============== STATISTICAL TESTING =============== STATISTICAL TESTING ==================

%% PCA All Data:
[testResults{iComb,1}.pca.coef,~,~,~,testResults{iComb,1}.pca.explain] = ...
    pca(testResults{iComb,1}.normAggMat(:,1:end));

%% MVRG All Data:

if smactinFlag.synthData(iMov) == 1 % Note: synthData == 1 case is obsolete. TN20191105
    
    try
        [testResults{iComb,1}.mvrgs.coef,testResults{iComb,1}.mvrgs.covarResponse, ...
            testResults{iComb,1}.mvrgs.residuals,testResults{iComb,1}.mvrgs.covarParam] = ...
            mvregress(aggMat2Analyze(:,1 : numSpklProp), ...
            aggMat2Analyze(:,numSpklProp + 1 : end),'maxIter',100);
    catch
    end
    
elseif smactinFlag.synthData(iMov) == 0
    
    % In case we pull out sm without any speckles around it, MVRG
    % calculation between speckles and sm will be between NaN values
    % and sm properties (function will not run!) Thus, we need to check if
    % speckles value are all NaN. If yes (i.e. no speckle around sm),
    % MVRG calculation is skipped, NaN is displayed as a flag that MVRG
    % was not run. - Tra Ngo 2019/07/05
    % This is obsolete due to implementation of completeFlag, only complete
    % dataset (matrix w/o NaN) will be pushed through mvrg - TN 2020/12/13
    
    sumNormAggMat = sum(testResults{iComb, 1}.normAggMat,'omitnan');
    
    % in normAggMat:
    % colume 1:numIndVars contains information for speckles
    % colume numIndVars + 1: end - 1 contains information for sm
    
    % if any normAggMat values are only NaN, output testResults as NaN.
    if any(sumNormAggMat==0) % sumNormAggMat(1,iColsumNormAggMat) == 0
        
        testResults{iComb,1}.mvrgs = mvrgTestNaN();
        
        % Report on whether the problem of NaN comes from FSM or SM:
        iColsumNormAggMat = find(sumNormAggMat == 0);
        if ismember(iColsumNormAggMat,1:numSpklProp)
            testResults{iComb,1}.mvrgs.nanFlag = 'NaN FSM. MVRG not run.';
        else
            testResults{iComb,1}.mvrgs.nanFlag = 'NaN SM. MVRG not run.';
        end
        
    else
        testResults{iComb,1} = fillStatTest(testResults{iComb,1}, propVect, ...
            numSpklProp, mvrgAggMatFlag, minNumData, randMvrgInfo);
        
    end % ( if any(sumNormAggMat==0) )
    
    
else
    
    error('Invalid input for handling synthetic data.')
    
end % ( if smactinFlag.synthData(iMov) == 1 )

%%  Statistical Test by Motion Types (AKA "Trajectory Classes")
%% MVRG each motion type:
testResults{iComb,2} = fillStatTest(testResults{iComb,2}, propVect, numSpklProp, mvrgAggMatFlag, minNumData, randMvrgInfo);

%% PCA each motion type:
for iDiff = 0 : 3
    [testResults{iComb,2}.pca{iDiff + 1}.coef,~,~,~,testResults{iComb,2}.pca{iDiff + 1}.explain] = ...
        pca(testResults{iComb, 2}.normAggMat{1, iDiff+1}(:,:));
    % 20191213: matlab auto normalization has different PCA results as current analysis
end

%% Split into 2 if needed
if dividingDatFlag
    testResults_actinCutOff{1,1} = fillStatTest(testResults_actinCutOff{1,1}, propVect, numSpklProp, mvrgAggMatFlag, minNumData, randMvrgInfo);
    testResults_actinCutOff{1,2} = fillStatTest(testResults_actinCutOff{1,2}, propVect, numSpklProp, mvrgAggMatFlag, minNumData, randMvrgInfo);
    testResults_actinCutOff{2,1} = fillStatTest(testResults_actinCutOff{2,1}, propVect, numSpklProp, mvrgAggMatFlag, minNumData, randMvrgInfo);
    testResults_actinCutOff{2,2} = fillStatTest(testResults_actinCutOff{2,2}, propVect, numSpklProp, mvrgAggMatFlag, minNumData, randMvrgInfo);
end

testResults_actinCutOff = testResults_actinCutOff';
testResults = testResults';
end


%% Sub-function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mvrgTestResult] = mvrgTestSorter(data1, data2, minNumData)
%MVRGTESTSORTER        Data analysis of multiple single-molecule and actin movies
%
%SYNOPSIS: [mvrgTestResult] = mvrgTestSorter(data, testResult)
%
%INPUT:        data: data matrix of single-molecule and actin properties in
%                    either normAggMat or aggMat format
%
%        testResult: cell that will contains MVRG results:
%                    .coef   .covarResponse   .residuals   .covarParam
%
%OUTPUT: mvrgTestResult: cell that contains MVRG results:
%                    .coef   .covarResponse   .residuals   .covarParam
%Tra Ngo, May 2020

% Initiate
mvrgTestResult.coef = [];
mvrgTestResult.covarResponse = [];
mvrgTestResult.residuals = [];
mvrgTestResult.covarParam = [];

[data1, nanColIndx1] = removeNanCol(data1); %design matrix (X, actin)
[data2, nanColIndx2] = removeNanCol(data2); %response matrix (Y, SM)

%% MVRG
if isempty(data1) || isempty(data2)  || min(size(data1,1),size(data2,1)) < minNumData
    [mvrgTestResult] = mvrgTestNaN();
else
    
    try
        [mvrgTestResult.coef, mvrgTestResult.covarResponse, ...
            mvrgTestResult.residuals, mvrgTestResult.covarParam] = ...
            mvregress(data1, data2);
    catch ME
        [mvrgTestResult] = ME;
    end
    
end

%% Fill coefficients row with NaN if input of design matrix (speckle properties) is NaN
if isstruct(mvrgTestResult) && ~isfield(mvrgTestResult,'nanFlag') && ~isempty(nanColIndx1) && (~isempty(data1) || ~isempty(data2))
    
    for i = 1:length(nanColIndx1)
        % if there is a col k of X (speckle, data1) that is NaN,
        % insert NaN in row k of coef matrix
        k = nanColIndx1(i);
        insertRow = nan(1,size(mvrgTestResult.coef,2));
        mvrgTestResult.coef = ...
            [mvrgTestResult.coef(1:(k-1),:); insertRow; mvrgTestResult.coef(k:end,:)];
        
        [nRowCoef,nColCoef] = size(mvrgTestResult.coef);
        
        % Insert NaN into covarParam
        for iIns = 1:nColCoef
            insPost = nRowCoef*iIns - 1;
            insertCol = nan(size(mvrgTestResult.covarParam,1),1);
            mvrgTestResult.covarParam = ...
                [mvrgTestResult.covarParam(:,1:(insPost-1)), insertCol, ...
                mvrgTestResult.covarParam(:,insPost:end)];
            insertRow = nan(1,size(mvrgTestResult.covarParam,2));
            mvrgTestResult.covarParam = ...
                [mvrgTestResult.covarParam(1:(insPost-1),:); insertRow; ...
                mvrgTestResult.covarParam(insPost:end,:)];
        end
    end
end
%% Fill coefficients col with NaN if input of response matrix (SM properties) is NaN
if isstruct(mvrgTestResult) && ~isfield(mvrgTestResult,'nanFlag') && ~isempty(nanColIndx2) && (~isempty(data1) || ~isempty(data2))
    % if there is a col k that is NaN, insert NaN in col k of coef matrix
    for i = 1:length(nanColIndx2)
        k = nanColIndx2(i);
        insertCol = nan(size(mvrgTestResult.coef,1),1);
        mvrgTestResult.coef = ...
            [mvrgTestResult.coef(:,1:(k-1)), insertCol, mvrgTestResult.coef(:,k:end)];
        clear insertCol
        nRowCoef = size(mvrgTestResult.coef,1);
        
        % Insert NaN into covarResponse
        insertCol = nan(size(mvrgTestResult.covarResponse,2),1);
        mvrgTestResult.covarResponse = ...
            [mvrgTestResult.covarResponse(:,1:(k-1)), insertCol, ...
            mvrgTestResult.covarResponse(:,k:end)];
        insertRow = nan(1,size(mvrgTestResult.covarResponse,2));
        mvrgTestResult.covarResponse = ...
            [mvrgTestResult.covarResponse(1:(k-1),:); insertRow; ...
            mvrgTestResult.covarResponse(k:end,:)];
        
        % Insert NaN into residuals
        insertCol = nan(size(mvrgTestResult.residuals,1),1);
        mvrgTestResult.residuals = ...
            [mvrgTestResult.residuals(:,1:(k-1)), insertCol, ...
            mvrgTestResult.residuals(:,k:end)];
        
        % Insert NaN into covarParam
        insertPosition = (k-1)*nRowCoef;
        insertCol = nan(size(mvrgTestResult.covarParam,2), nRowCoef);
        mvrgTestResult.covarParam = ...
            [mvrgTestResult.covarParam(:,1:insertPosition), insertCol, ...
            mvrgTestResult.covarParam(:,insertPosition+1:end)];
        insertRow = nan(nRowCoef,size(mvrgTestResult.covarParam,2));
        mvrgTestResult.covarParam = ...
            [mvrgTestResult.covarParam(1:insertPosition,:); insertRow; ...
            mvrgTestResult.covarParam(insertPosition+1:end,:)];
    end
end

end

function [dataOut, nanColIndx] = removeNanCol(dataIn)
%REMOVENANCOL       check which column inside a matrix is composed of NaN entirely and removew that column

if isempty(dataIn)
    nanColIndx = [];
else
    nanColIndx = find(~sum(~isnan(dataIn),[1])); % find which column composed of just NaN values entirely
    dataIn(:,nanColIndx) = []; % remove the column of NaN
end

dataOut = dataIn;
end

function [mvrgTestResult] = mvrgTestNaN()
%MVRGTESTNAN        Fill data analysis result of multiple single-molecule and actin movies with NaN
% if there is no information for either sm or speckles, mvrg
% calculation will not run!
%SYNOPSIS: [mvrgTestResult] = mvrgTestSorter(data, testResult)
%
%INPUT:        data: data matrix of single-molecule and actin properties in
%                    either normAggMat or aggMat format
%
%        testResult: cell that will contains MVRG results:
%                    .coef   .covarResponse   .residuals   .covarParam
%
%OUTPUT: mvrgTestResult: cell that contains MVRG results:
%                    .coef   .covarResponse   .residuals   .covarParam
%Tra Ngo, May 2020
mvrgTestResult.coef = NaN;
mvrgTestResult.covarResponse = NaN;
mvrgTestResult.residuals = NaN;
mvrgTestResult.covarParam = NaN;
mvrgTestResult.nanFlag = 'MVRG not ran. Either (1) SM or FSM (or their normalized data) is NaN or (2) not enough data by min requirement.';
end

function [normDataMat] = normalizeMat(dataMat,methodNum)
%NORMALIZEAGGMAT    takes a matrix (.aggMat) and calculate normalized data(.normAggMat)
%INPUT: dataMat (matrix) : col = properties; row = observation
%       methodNum (integer): 1 or 2
%
%OUTPUT: normDataMat (matrix) : col = properties; row = observation
%
%Tra Ngo, July 2020
if ~exist('methodNum','var') || isempty(methodNum)
    methodNum = 2;
end
[nObservation, nProP] = size(dataMat);

switch methodNum
    case 1
        %% Method 1: somehow gain more NaN in results than method 2
        oneVect = ones(nObservation,1);
        meanVect = nanmean(dataMat);
        stdVect = nanstd(dataMat); stdDiagMat = diag(stdVect);
        
        normDataMat = (dataMat - oneVect*meanVect)/(stdDiagMat); % Standard normalization of dataMat
    case 2
        %% Method 2: default
        normMatAll = [];
        for iProp = 1 : nProP
            normProp = (dataMat(:,iProp) - nanmean(dataMat(:,iProp))) / ...
                nanstd(dataMat(:,iProp));
            normMatAll = horzcat(normMatAll,normProp);
        end
        normDataMat = normMatAll;
    otherwise
        disp('Method undefined')
end
end

function [testResultsOut] = fillStatTest(testResults, propVect, numSpklProp, mvrgAggMatFlag, minNumData, randMvrgInfo)
%FILLSTATTEST: (sub-function) input structure with .aggMat and .normAggMat;
%INPUT: testResult without fields .mvrgs and .mvrgRand
%       propVect: vector of properties of interest
%       numSpklProp: number of speckle properties
%       mvrgAggMatFlag: flag indicating if normalized data should be used
%       minNumData: (integer) minimum number of data point that should
%                   exist to run MVRG.
%       randMvrgInfo: (1 x 2 cell) 1st cell is logic flag indicating
%                      whether randomization should be performed. 2nd cell
%                      is the number of randomization should be performed.
%OUTPUT: structure with .mvrgs and .mvrgRand (optional) added.
testResultsOut = testResults;
randMvrgFlag = randMvrgInfo{1};
numRand = randMvrgInfo{2};

if isa(testResultsOut.aggMat,'cell')
    nType = 4; % if .aggMat is cell containing 4 matrix for 4 motion types
else
    nType = 1;
end

for iDiff = 1 : nType
    
    % Calculate X
    
    if mvrgAggMatFlag
        
        try 
            compltIndx = testResultsOut.completeFlag{iDiff};
            compltData = testResultsOut.aggMat{iDiff}(compltIndx,propVect);
            inlierIndx = testResultsOut.outlierFlag{iDiff};
            dataCurr = compltData(inlierIndx,:);

        catch
            compltIndx = testResultsOut.completeFlag;
            compltData = testResultsOut.aggMat(compltIndx,propVect);
            inlierIndx = testResultsOut.outlierFlag;
            dataCurr = compltData(inlierIndx,:);
        end
        
        lengthX = size(dataCurr,1);
        oneVect = ones(lengthX,1);
        X = [oneVect dataCurr(:,1 : numSpklProp)];
    else
        
        try
            dataCurr = testResultsOut.normAggMat{iDiff};
        catch
            dataCurr = testResultsOut.normAggMat;
        end
        
        X = dataCurr(:,1 : numSpklProp);
        
    end % (if mvrgAggMatFlag)
    
    % Calculate Y
    
    Y = dataCurr(:,numSpklProp + 1 : end);
    
    % MVRG
    
    if nType == 4
        testResultsOut.mvrgs{iDiff} = mvrgTestSorter(X,Y,minNumData);
    elseif nType == 1
        testResultsOut.mvrgs = mvrgTestSorter(X,Y,minNumData);
    else
        error("Case not consided in fillStatTest!");
    end
    
    
    % MVRG of randomized
    if randMvrgFlag
        
        lenX = size(X,1); % number of observation in X
        
        for iRand = 1:numRand
            % Randomize X (speckle)
            randVect = randperm(lenX,lenX); % create randomized index vector
            randX = X(randVect,:);
            
            if nType == 4
                testResultsOut.mvrgRand{iDiff}{iRand,1} = mvrgTestSorter(randX,Y,minNumData);
            elseif nType == 1
                testResultsOut.mvrgRand{iRand,1} = mvrgTestSorter(randX,Y,minNumData);
            else
                error("Case not consided in fillStatTest!");
            end
        end % (for iRand = 1:100)
        
    end % (if randMvrgFlag)
    
end % (for iDiff = 1 : nType)

end

function [dataMatOut, completeFlag, normAggMatInlier, outlierLabelVect, normDataComplete] = ...
    processCompltOutliersFrMat(dataMat,propVect, neighborNum, iterThres, maxIter, k, indxColRmv, numCol)
% PROCESSCOMPLTOUTLIERSFRMAT: process dataMat (aggMat) to give label of complete data point and label of outlier (outlier detection is not needed here anymore as of 20210722)
% INPUT:  dataMat (matrix) of  all properties (cols) and all observations (rows)
%
%         propVect (vector) of cols in dataMat that are of interests
%
%         neighborNum (interger) number of neighbor to calculate outliers
%
%         iterThres (optional number) .default = 0.1 (10%) ratio of
%                       tolerance to stop iterating outlier detection.
%
%         maxIter (optional integer) .default = 4, maximum number of iteration
%                       for outlier detection.
%
%           k : (optional, integer) (default: [], so that value is 
%               determined by labelOutliersFrMat.m's default input).
%                   k= -1: outlier detection will not be run
%                   k > 0: roughly, observations that are k*sigma away from 
%                          the mean  will be considered outliers.
%
%          indxColRmv: (obsolete) indexing which col is SM diffusion
%                     coefficient, because in outlier detection, we do not want 
%                     to count identical diffusion coefficients as neighbors.
%
%           numCol: (integer) indicating the number of column in an
%                   expected dataMat. Used in case dataMat is empty and we need 
%                   to initialize outputs in the correct 0x(numCol) format.
%                     i.e.: numCol = 14; empty dataMat is set to 0x14
%
% OUTPUT: dataMatOut is equivalent to aggMat
%
%         normAggMatInlier is  equivalent to normAggMat of complete & inlier-only aggMat
%
% assumption of dataMat: 3rd col is speckle density (0 for no matching)
% this function looks at this column, if there's no matching, matchFlag
% output NaN for that sm observation.
%% Outputs
completeFlag = []; %zeros(0,1);
dataMatOut = dataMat;
normAggMatInlier = []; %zeros(0,2);
outlierLabelVect = []; %zeros(0,1);
normDataComplete = []; %zeros(0,2);
%% Inputs control & constant
if isempty(dataMat)
    % warning('Input dataMat is empty.');
    dataMat = zeros(0,numCol);
    dataMatOut = dataMat;
    %return
end
   
if ~exist("iterThres","var") || isempty(iterThres)
    iterThres = 0.1; % 10%
end
if ~exist("maxIter","var") || isempty(maxIter)
    maxIter = 4;
end
if ~exist('k','var')
    k = [];
end
colSpklDens = 3;

%% Match Flag: & replace speckle density of non-matching with NaN
%1 for observations matching sm and at least 1 speckle
%0 for observations matching sm and NO speckle
% if size(dataMat,2) < colSpklDens
%     warning('dataMat (aggMat) do not contain info for speckle density')
% elseif isnan(mean(dataMat(:,colSpklDens)))
%     warning('Speckle density contains NaN. Results matrix do not search for any actin density of 0 left.')
% else
%     matchFlag = dataMat(:,colSpklDens) ~= 0;
%     dataMatOut(~matchFlag,colSpklDens) = NaN;
% end

if size(dataMat,2) < colSpklDens
    warning('dataMat (aggMat) do not contain info for speckle density')
else
    matchFlag = dataMat(:,colSpklDens) == 0;
    dataMatOut(matchFlag,colSpklDens) = NaN;
end


%% Complete Flag: from aggMat contains properties of interest
%1 = observations having NO NaN; 0 = observations having any NaN
matOfInterest = dataMatOut(:,propVect);
completeFlag = ~isnan(mean(matOfInterest,2));
completeDatOfInterest = dataMatOut(completeFlag, propVect); % AggMat with incomplete data removed

%% Normalized aggregate matrix of matched properties of nn combinations
normDataComplete = normalizeMat(completeDatOfInterest,2);

%% Outliers: Test if detection is good by removing outliers -> reapply
if k ~= -1
    % this process is iterated until (# of outliers removed in rorund n+1)/(# of
    % outliers removed in round n) < 10% (or defined by iterThres); maximum
    % number of iteraction for outlier detection is defined by maxIter
    iTer = 1;
    % Detect outliers (0s) among the normalized, matched sm-spkl observations
    outlierLabelVectTemp1 = labelOutliersFrMat(normDataComplete, neighborNum, k, indxColRmv);
    % % indexInlier = find(outlierLabelVectTemp1);
    indexInlier = 1:size(normDataComplete,1);
    normDatComplInlTemp1 = normDataComplete(outlierLabelVectTemp1,:); % remove outliers from normalized data
    nOutPrevRound = sum(~outlierLabelVectTemp1); % if no outlier is detected, no subsequent iteration is done
    nOutPostRound = inf;
    convRat = nOutPostRound/nOutPrevRound;
    
    normDatComplInlTemp2 = normDatComplInlTemp1;
    outlierLabelVectTemp2 = outlierLabelVectTemp1;
    
    while nOutPrevRound ~= 0 && convRat >= iterThres && iTer <= maxIter
        % Copy result from previous iteration & report index of inliers
        outlierLabelVectPrev = outlierLabelVectTemp2;
        normDatComplInlTempPrev = normDatComplInlTemp2;
        indexInlier = indexInlier(outlierLabelVectPrev);
        convRatPrev = convRat;
        
        % Next iteration of outlier detection on normalize data (with outliers removed)
        outlierLabelVectTemp2 = labelOutliersFrMat(normDatComplInlTempPrev, neighborNum, k, indxColRmv);
        normDatComplInlTemp2 = normDatComplInlTempPrev(outlierLabelVectTemp2 ,:);
        
        % Calculate ratio of number outliers (to cal. tolerance for res convergence)
        nOutPostRound = sum(~outlierLabelVectTemp2);
        convRat = nOutPostRound/nOutPrevRound;
        
        % Copy result for next iteration
        nOutPrevRound = nOutPostRound;
        
        iTer = 1+iTer;
        %{
         % Test if outliers detection is good VISUALLY
%         x = normDataComplete(:,3); % spkl
%         y = normDataComplete(:,4); % sm
%         gscatter(x,y,outlierLabelVect,'br','xo'); % 'br' indicates blue-red
%         xlabel('spkl lifetime');
%         ylabel('sm mean amp');
        %}
        
        % Save outlier label
    end
    
    disp(['outlier detection completed at iteration: ' num2str(iTer-1) ' with convergence ratio (if continue with next iteration) = ' num2str(convRat)]); % this display is a little buggy (if not enter while-loop, iTer = 0)
    % if function for outlier detection is ran once and it already provides
    % enough outliers, num2str(iTer-1) = 0 (do not enter for-loop AKA no outliers)
    % or 1 (enter for-loop once)
    outlierLabelVect = false(size(normDataComplete,1),1);
    outlierLabelVect(indexInlier) = true;
else
    outlierLabelVect = true(size(normDataComplete,1),1);
end

%% Remove the outliers and re-normalize the data
aggMatCompleteInlier = completeDatOfInterest(outlierLabelVect,:);
normAggMatInlier = normalizeMat(aggMatCompleteInlier,2);
end

function [abvAggMat, abvPos, belAggMat, belPos] = getCompltInlierIndex(dataMat, arrComplete, arrInlier, arrPos, actinCutOff, spklDensityCol, numCol)
%GETCOMPLTINLIERINDEX: take testResults structure and return index of "above" or "below" for the field 'aggMat'
%INPUT:   dataMat       : (matrix) data, col = properties, row = observations 
%                         i.e.: testResults{1,1}.aggMat
%         arrComplete   : vector of logic indicating "1" = complete data, "0"
%                         = incomplete data 
%                         i.e.: testResults{1,1}.completeFlag 
%         arrInlier     : vector of logic indicating "1" = inliers, "0" = outliers
%                         i.e.: testResults{1,1}.outlierFlag
%         arrPos        : nx2 vector  of mean position of sm for each datapoint
%                         i.e.: testResults{1,1}.smMeanPos
%         actinCutOff   : i.e.: 0.01
%         spklDensityCol: integer indicates the colume in aggMat that
%                         contains info of property that data is split by
%                         (i.e. speckle density = 3)
if isempty(dataMat)
    disp('dataMat is empty')
    dataMat = zeros(0,numCol);
end
dataComplt = dataMat(arrComplete,:);
dataCompltInl = dataComplt(arrInlier,:);

indxAbv = dataCompltInl(:,spklDensityCol) > actinCutOff;
indxBel = ~indxAbv;

abvAggMat = dataCompltInl(indxAbv,:);
belAggMat = dataCompltInl(indxBel,:);

% note for if which column is "above cutoff" ever become ambiguous:
% testResults_actinCutOff{1,1}.note = 'Above';
try
    abvPos = arrPos(indxAbv,:);
    belPos = arrPos(indxBel,:);
catch
    abvPos = zeros(0,2);
    belPos = zeros(0,2);
end

end