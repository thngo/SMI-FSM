function [inputBatchName, pSigIndCellArr, outlierIndxArr, sigRes, mvrgCoefStructArr] = multiCompareIndCellsMVRG(spekLabelBase, smLabelBase, displSpacing, varargin)
%MULTICOMPAREINDCELLSMVRG: takes different testResult containing coefficients of multivariate regression of individual cells from different population (i.e. treatment conditions) and return significance comparisons.
%
%SYNOPSIS: [inputBatchName, pSigIndCellArr, outlierIndxArr, sigRes, mvrgCoefStructArr] = multiCompareIndCellsMVRG(spekLabelBase, smLabelBase, varargin)
%
%INPUT   :  spekLabelBase: (cell array) names of speckle properties of interest.
%               i.e.: spekLabelBase = {'Speckle Speed', 'Speckle Local
%               Density','Speckle Lifetime'};
%
%           smLabelBase  : (cell array) names of SM properties of interest.
%               i.e.: smLabelBase = { 'SM Mean Amplitude','SM Diff Coef',
%               'SM Local Density','smDiffusionRadius', 'SM Aggreg State'};
%
%           displSpacing : allow user to input spacing parameters for final  
%                          output notBoxPlot figure [spaceBtwEachBoxPlot, spaceBtwEachGroup]. 
%
%           vargarin     : Each input is a cell array of individual cells
%                          testResult structure (which contain MVRG coefs).
%                          If first cell array (3rd function input) is
%                          logical value "true", then function expects only
%                          2 treatment populations (2 testResults) and
%                          perform ttest to compare, else, anovann is used.
%                          Default comparison test is ANOVA. 
%
% Example for 3 and more (run ANOVA test):
% [inputBatchName, pSigIndCellArr, outlierIndxArr, sigRes] = multiCompareIndCellsMVRG(spekLabelBase, smLabelBase, displSpacing, ...
%     strInl_1.testResultSplitAbv, strInl_2.testResultSplitAbv, strInl_3.testResultSplitAbv, strInl_4.testResultSplitAbv);
%
% Example for 2 (run ttest2) -- 3rd input MUST be true.
%[inputBatchName, pSigIndCellArr, outlierIndxArr, sigRes] = multiCompareIndCellsMVRG(spekLabelBase, smLabelBase, displSpacing, ...
%    true, strInl_1.testResultSplitAbv, strInl_2.testResultSplitAbv);
%
%OUTPUT  : function outputs folders containing different notBoxPlot figures
%          of MVRG coefficients.
%               For individual populations MVRG coefficients visualization:
%                  folders named "indCellsMVRG#_(name of input variable)"
%               For comparing between the different populations's MVRG
%                  coefficients: 1 folder named "compareSigAnova" or
%                  "compareSigTtest".
%
%          inputBatchName: array of cells, each containing string of
%                          variable names of different treatment population
%                          input. (i.e. "testResultsInd_C1_aboveTh")
%                          Note that if an input is direct access (i.e.
%                          a{1}), this string will be empty.
%
%          pSigIndCellArr: array of cells containing pSig of each population.
%
%          outlierIndxArr: array of cells containing outliers of each population.
%
%          sigRes        : either output of multiCompareAnovaIndCellsMVRG.m if AnovaN;
%                          or structure contains 4 outputs of compareSigIndCellsMVRG.m if t-test.
%
%          mvrgCoefStructArr: (n x 1 cell array) each cell contains MVRG
%                          coefficients from each batch of data that was input.
%                          Each cell contain a structure with the following fields:
%                                  .Immobile
%                                  .Confined
%                                  .Free
%                                  .Directed
%                                  .All
%                          Which contains fields corresponding to the
%                          dependpent variables from MVRG. Inside these
%                          field is a matrix for MVRG coefficients, with
%                          row = individual cell ROI (movie) and col =
%                          indpendent variables from MVRG, whose names are
%                          listed in the field "mvrgCoefStructArr.speckleLabel"
%
% Warning: Please always run this function in a new empty directory, because
% this function create folders with set names and output figures into those folder!
%
% Tra Ngo, April 2021
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


%% Input:
if nargin < 4
    error('Function expects at least 4 inputs, with the first 2 inputs being labels of speckle and single-molecule properties.')
end

if isempty(displSpacing)
    displSpacing = [0.4 5];
end

numBatch = nargin - 3; % number of testResult input
%% Output:
inputBatchName = cell(numBatch,1);
for i = 1:numBatch
    inputBatchName{i} = inputname(i+3);
end

%% Put variable inputs into cell array
batchInputArr = varargin';

%% Determine compare test: anovan? or ttest?
if islogical(batchInputArr{1}) && batchInputArr{1}
    if numBatch == 3
        ttestFlag = true;
        
        % remove logic flag from batchInputArr, excessive name from
        % inputBatchName, get the true number of treatment population being
        % tested (i.e. numBatch):
        
        batchInputArr(1) = [];
        inputBatchName(1) = [];
        numBatch = numBatch-1;
    else
        error('Unexpected number of treatment population if ttest is requested. Expecting only 2 treatment populations.')
    end
elseif islogical(batchInputArr{1}) && ~batchInputArr{1}
    error('Unexpected input: 3rd input cannot be logical value "false".')
else
    ttestFlag = false; % (default behavior)
end

%% Other outputs:
pSigIndCellArr = cell(numBatch,1);
outlierIndxArr = cell(numBatch,1);
mvrgCoefStructArr = cell(numBatch,1);

%% Check significance or individual population:

for iBch = 1: numBatch
    % batchInputArr = varargin{iBch};
    folderName = ['indCellsMVRG' num2str(iBch) '_' inputBatchName{iBch}];
    mkdir(folderName)
    cd(folderName)
    [pSigIndCellArr{iBch}, outlierIndxArr{iBch}, ~, mvrgCoefStructArr{iBch}] = ...
        explainIndCellsMVRG(batchInputArr{iBch}, spekLabelBase, smLabelBase);
    cd ..
end


%% Compare between populations:
if numBatch ~= 1
    switch ttestFlag
        
        case false % if ttestFlag is false, then perform compare using anovan
            
            mkdir('compareSigAnova')
            cd('compareSigAnova')
            [sigRes] = multiCompareAnovaIndCellsMVRG(batchInputArr,...
                pSigIndCellArr, outlierIndxArr, spekLabelBase, smLabelBase, displSpacing);
            cd ..
            
        case true % if ttestFlag is true, then perform compare using ttest
            mkdir('compareSigTtest')
            cd('compareSigTtest')
            [pCompareTtestR,pCompareTtestL,outlierIndx1,outlierIndx2] = ...
                compareSigIndCellsMVRG(batchInputArr{1},batchInputArr{2},...
                pSigIndCellArr{1},pSigIndCellArr{2},...
                spekLabelBase, smLabelBase);
            cd ..
            
            sigRes.pCompareTtestR = pCompareTtestR;
            sigRes.pCompareTtestL = pCompareTtestL;
            sigRes.outlierIndx1 = outlierIndx1;
            sigRes.outlierIndx2 = outlierIndx2;
            
    end % (switch ttestFlag)
else
    sigRes = 'No comparison. Only 1 population was input.';
end % (if numBatch == 1) then don't do any compare

%% If something comes out as different, implement multCompare to find what is different:

end