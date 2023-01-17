function outputStruct = runSmactinAggMovExplain(tassm,matchindx,smactinFlag, propSpk, propSm, ...
    maxNumPropVect, mvrgAggMatFlag, cutOffInput, outlierParam, spekLabelBase, smLabelBase, sigThreshInput)
%RUNSMACTINAGGMOVEXPLAIN takes data matrices of SM - speckles pairs, perform MVRG on them, and report statistical significances of the regression coefficients.
%
%SYNOPSIS: outputStruct = runSmactinAggMovExplain(tassm,matchindx,smactinFlag, propSpk, propSm, ...
%                         maxNumPropVect, mvrgAggMatFlag, cutOffInput, outlierParam, spekLabelBase, ...
%                         smLabelBase, sigThreshInput)
%
%PROCESS: runSmactinAggMovExplain.m takes individual data matrices from tassm, 
% call smactinAggMov.m for individual cells, split its result into "above" and 
% "below" populations, and run subsequent explaination functions of the MVRGs.
%
%INPUT: tassm           : (n x 1 cell array of cells) Each cell contains 3 x 1 
%                         cells, containing data matrix of SM tracklets and 
%                         speckle properties:
%                       Cell {1,1}: (k observations x 9 properties) Data 
%                       matrix of speckle properties with columns as:
%                               1: Local speckle displacement magnitude
%                               2: Speckle co-movement
%                               3: Local speckle density 
%                               4: Speckle lifetime
%                               5: Local speckle intensity
%                               6: Local speckle background intensity
%                               7: speckle corrected intensity
%                               8: speckle corrected background intensity
%                               9: Local speckle average displacement magnitude
%                       Cell {2,1}: (k observations x 7 properties) Data 
%                       matrix of SM properties with columns as:
%                               1: SM net speed
%                               2: SM intensity
%                               3: SM diffusion coefficient
%                               4: SM density
%                               5: SM spatial span
%                               6: SM apparent assembly state (oligomerization state)
%                               7: SM diffusion type
%                       See smactinMat.m for details.
%
%        matchindx      : indices of speckles,ks, or sms matching to an object;
%                         empty if many neighbors matching is performed.
%                       See smactinMat.m for details.
%
%        smactinFlag    : (structure) flags describing the method by which 
%                         SM tracklets and actin speckles are associated
%                         spatially; contains the following field: 
%                               .combFlag: integer with value either 7 (to 
%                                 match SM to speckles based on average SM
%                                 position) or 10 (to match SM to speckles
%                                 based on frame-by-frame SM position).
%                               .randFlag: (integer) must put value > 0.
%
%        propSpk        : (vector of interger) of speckle properties in tassm
%                         collected for normAggMat. (optional) If input as
%                         empty [] then default to [1 3 4].
%                         See smactinAggMov.m for details.
%
%        propSm         : (vector of interger) of SM properties in tassm
%                         collected for normAggMat. (optional) If input as
%                         empty [] then default to [2 3 4].
%                         See smactinAggMov.m for details.
%
%        maxNumPropVect : (1x2 vector) number of speckle property in tassm
%                         and number of SM property in tassm. 
%
%        mvrgAggMatFlag : (logic, optional) flag indicating whether
%                         normalization should be performed before MVRG.
%                         default. 0 = MVRG on normAggMat. (normalized data)
%                                  1 = MVRG on aggMat. (original data;
%                                      cutting through intercept [oneVect X-data])
%
%        cutOffInput    : (1x2 vector, optional) threshold at which data is 
%                        split into 2 subsets before any MVRG and
%                        statistical analysis.
%                        (default) empty [] -> data not split into 2.
%                        1st entry: threshold dividing the data into 2.
%                        2nd entry: the column in .aggMat that the threshold 
%                        is scanned within, to split the data into 2.
%                        Example: cutOffInput = [0.01 3];
%
%        outlierParam   : (4 x 1 vector, optional) Parameters for outlier
%                        detection process, containing in order:
%                        [neighborNum, iterThres, maxIter, k]. If [] is
%                        input, default to [2 0.1 4 -1] (empirically determined).
%                        See smactinAggMov.m for details.
%
%        spekLabelBase  : (1 x n cell array) each cell contains name of 
%                        speckle property to be plotted on the x-axis.
%                        See explainIndCellsMvrgWithRand.m for details.
%
%        smLabelBase    : (1 x n cell array) each cell contains name of SM 
%                        property for plot title.
%                        See explainIndCellsMvrgWithRand.m for details.
%
%        sigThreshInput : (1x2 double vector, optional) alpha thresholds 
%                        that determine whether p-values are significant. 
%                        [alpha-value for ttest , alpha-value for randomization test] 
%                        Optional. If empty [] was input, default to [0.05, 0.05].
%
%
% ------- example inputs -------
% smactinFlag = struct('combFlag',10,'randFlag',100, 'match',2, 'synthData',0);
% propSpk = [1 3 4 5]; % speed, density, lifetime, ilmax
% propSm = [2 3 4]; % mean amp, diff coef, density
% maxNumPropVect = [9 6]; % 9 speckle props, 6 sm props
% mvrgAggMatFlag = 0; % MVRG done on normAggMat
% cutOffInput = [0.01 3]; % threshold on actin density at 0.01
% outlierParam = [2 0.1 4 -1]; % -1 for not running any outlier detection
% sigThreshInput = [0.05 0.05];
% spekLabelBase = {'Speckle Speed','Speckle Local Density','Speckle Lifetime'};
% smLabelBase = {'SM Mean Amplitude','SM Diff Coef','SM Local Density'};
% ------------------------------
%
%OUPUT: Function outputs 2 folders: "above" and "below" that contain
%       box-plot of regression coefficients from MVRG performed on data 
%       above and below the user-input threshold, respectively.
%
%       outputStruct: (structure) results of MVRG and relevant statistical 
%                     analyses relating to comparing MVRG coefficients;
%                     contains the following fields:
%                     .testResultAndSplit: MVRG results of all input data
%                               and of any subsets of data below and above 
%                               the user-input threshold.
%                               See smactinAggMov.m for details.
%                     .trajClassAndInputParam: {1,1} contains SM track's
%                               diffusion type (See smactinAggMov.m for
%                               details); {1,2} contains input parameters
%                               of this function.
%                               See smactinAggMov.m for details.
%                     .testResultSplitAbv: MVRG results of subset of data
%                               above the user-input threshold.
%                               See smactinAggMov.m for details.
%                     .testResultSplitBel: MVRG results of subset of data
%                               below the user-input threshold.
%                               See smactinAggMov.m for details.
%                     .explainIndCellsMvrgWithRand_abv: significances for
%                               MVRG results of randomized data of
%                               individual cells, for data above the
%                               user-input threshold.
%                               See explainIndCellsMvrgWithRand.m for details.
%                     .explainIndCellsMvrgWithRand_bel: significances for
%                               MVRG results of randomized data of
%                               individual cells, for data below the
%                               user-input threshold.
%                               See explainIndCellsMvrgWithRand.m for details.
%                     .explainIndCellsMVRG_abv: significances and outliers
%                               of MVRG coefficients, for data above the
%                               user-input threshold.
%                               See explainIndCellsMVRG.m for details.
%                     .explainIndCellsMVRG_bel: significances and outliers
%                               of MVRG coefficients, for data below the
%                               user-input threshold.
%                               See explainIndCellsMVRG.m for details.
%
%GLOSSARY:  SMI   : single molecule imaging
%           SM    : single molecule
%           FSM   : fluorescent speckle microscopy
%           tassm : total attributes speckles and single molecules
%           MVRG  : multivariate regression
%
%Tra Ngo, Aug 2021
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

%% Inputs & initialize variables
len = length(tassm);
if isempty(matchindx)
    matchindx = cell(len,1);
end

testResult = cell(len,1);
trajClass = cell(len,1);
testResultSplit = cell(len,1);
inputParam = cell(len,1);

%% Run MVRG analysis for each cell ROI/movie
for i = 1:len
    [testResult{i,1},trajClass{i,1}, testResultSplit{i,1}, inputParam{i,1}] = ...
        smactinAggMov(tassm(i),matchindx(i),smactinFlag, ...
        propSpk, propSm, maxNumPropVect, mvrgAggMatFlag, cutOffInput, outlierParam); 
end

testResultSplitAbv = cell(len,1);
testResultSplitBel = cell(len,1);
for i = 1:len
    testResultSplitAbv{i,1} = testResultSplit{i,1}(:,1);
    testResultSplitBel{i,1} = testResultSplit{i,1}(:,2);
end


outputStruct.testResultAndSplit{1} = testResult;
outputStruct.testResultAndSplit{2} = testResultSplit;

outputStruct.trajClassAndInputParam{1} = trajClass;
outputStruct.trajClassAndInputParam{2} = inputParam;

outputStruct.testResultSplitAbv = testResultSplitAbv;
outputStruct.testResultSplitBel = testResultSplitBel;

[belowThreshFlag3D_abv, fraction3D_abv, avgPval_abv, coefMat4D_abv] = explainIndCellsMvrgWithRand(testResultSplitAbv, spekLabelBase, smLabelBase, sigThreshInput);
[belowThreshFlag3D_bel, fraction3D_bel, avgPval_bel, coefMat4D_bel] = explainIndCellsMvrgWithRand(testResultSplitBel, spekLabelBase, smLabelBase, sigThreshInput);

outputStruct.explainIndCellsMvrgWithRand_abv{1} = belowThreshFlag3D_abv; 
outputStruct.explainIndCellsMvrgWithRand_abv{2} = fraction3D_abv; 
outputStruct.explainIndCellsMvrgWithRand_abv{3} = avgPval_abv; 
outputStruct.explainIndCellsMvrgWithRand_abv{4} = coefMat4D_abv; 

outputStruct.explainIndCellsMvrgWithRand_bel{1} = belowThreshFlag3D_bel; 
outputStruct.explainIndCellsMvrgWithRand_bel{2} = fraction3D_bel; 
outputStruct.explainIndCellsMvrgWithRand_bel{3} = avgPval_bel; 
outputStruct.explainIndCellsMvrgWithRand_bel{4} = coefMat4D_bel; 

mkdir above; cd above
[pSigTtest_abv, outlierIndx_abv, randResult_abv, mvrgCoef_abv] = explainIndCellsMVRG(testResultSplitAbv, spekLabelBase, smLabelBase, sigThreshInput);
cd ..; mkdir below; cd below
[pSigTtest_bel, outlierIndx_bel, randResult_bel, mvrgCoef_bel] = explainIndCellsMVRG(testResultSplitBel, spekLabelBase, smLabelBase, sigThreshInput);
cd ..

outputStruct.explainIndCellsMVRG_abv{1} = pSigTtest_abv; 
outputStruct.explainIndCellsMVRG_abv{2} = outlierIndx_abv; 
outputStruct.explainIndCellsMVRG_abv{3} = randResult_abv; 
outputStruct.explainIndCellsMVRG_abv{4} = mvrgCoef_abv;

outputStruct.explainIndCellsMVRG_bel{1} = pSigTtest_bel; 
outputStruct.explainIndCellsMVRG_bel{2} = outlierIndx_bel; 
outputStruct.explainIndCellsMVRG_bel{3} = randResult_bel; 
outputStruct.explainIndCellsMVRG_bel{4} = mvrgCoef_bel;
end