function [belowThreshFlag3D, fractionFlagMat3D, avgPval, coefMat4D] = explainIndCellsMvrgWithRand(testResults, spekLabelBase, smLabelBase, sigThreshInput)
%EXPLAININDCELLSMVRGWITHRAND calculate significances for multivariate regression results of randomized data of individual cells
%
%SYNOPSIS:  [belowThreshFlag3D, fractionFlagMat3D, avgPval] = explainIndCellsMvrgWithRand(testResults, spekLabelBase, smLabelBase, sigThreshInput)
%
%INPUT     testResults  : Cell array of individual cell analysis results.
%                         This function cares about the multivariate
%                         regression output within this cell array. See
%                         smactinAggMov for details.
%
%          spekLabelBase: (1 x n cell array) each cell contains name of 
%                         speckle property to be plotted on the x-axis.
%          smLabelBase  : (1 x n cell array) each cell contains name of SM 
%                         property for plot title.
%               Choose appropriate speckle and SM properties, which were
%               used for MVRG's testResults and input in an order
%               consistent with testResults. 
%                   Example: if you previously perform MVRG with speckles
%                   columns [1 3 4] and with SM columns [2 3 4 6], then the
%                   following would be your input:
%                       spekLabelBase = {'Speckle Speed', 'Speckle Local Density','Speckle Lifetime'};
%                       smLabelBase = { 'SM Mean Amplitude','SM Diff Coef', 'SM Local Density', 'SM Aggreg State'};
%                   Remark: testResults currently can have the follwing:
%                     propVect1: (vector of interger) of speckle properties in tassm
%                     collected for normAggMat in testResults
%                       1: speckleSpeed
%                       2: speckleMvmtCohere (speckle co-movement)
%                       3: speckleDensity
%                       4: speckleLifetime
%                       5: speckle intensity (ilmax)
%                       6: speckle background intensity (ibkg)
%                       7: ilmaxCorrected
%                       8: ibkgCorrected
%                       9: speckle average speed
%                     propVect2: (vector of interger) of SM properties in tassm
%                     collected for normAggMat in testResults
%                       1: smNetSpeed
%                       2: smMeanAmp
%                       3: smDiffusionCoefficients
%                       4: smDensity
%                       5: smDiffusionRadius
%                       6: smAggregState (oligomerization state)
%                     See smactinAggMov.m for details.
%
%           sigThresh   : (1x2 double vector, optional) alpha thresholds 
%                        that determine whether p-values are significant. 
%                        [alpha-value for ttest , alpha-value for randomization test] 
%                        Optional. If empty [] was input, default to [0.05, 0.05].
%
%OUTPUT  belowThreshFlag3D  : 3D array of flags indicating which
%                             coefficient have significance p-values.
%                             (when percentage of times randomized data is
%                             deemed as significant by alpha value
%                             sigThreshInput(1) is smaller than
%                             sigThreshInput(2))
%                               1 = equal or fewer times than sigThreshInput(2)
%                               0 = more times than sigThreshInput(2)
%
%        fractionFlagMat3D  : indicates fraction of randomization (0 to 1)
%                             contains significant coefficients by ttest
%                             with alpha-value = alphaTtest.  
%
%        avgPval            : 3D array of average p-values of randomization 
%                             cases which was determined to be significant. 
%
%        coefMat4D          : (3D array) p-value from ttest of each
%                             randomization cases (added TN20220331)
%                             Dim as follow: (iSpk, iSMProp, iTraj, iRand)
%
% Dimensions are (number of speckle properties) by (number of SM
% properties) by (number of trajectory types). Trajectory types are:
% immobile, confined, free, directed and ALL.
%
%GLOSSARY:  SMI   : single molecule imaging
%           SM    : single molecule
%           FSM   : fluorescent speckle microscopy
%           tassm : total attributes speckles and single molecules
%           MVRG  : multivariate regression
%
%Tra Ngo, 08/2021
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
if nargin < 3 %input error
    error('Expected number of input is at least 3.');
end

if ~iscell(testResults{1})
    error('Input must contain more than 1 movie.');
end

numCell = length(testResults); % number of cells
numSmLabel = length(smLabelBase);
groupsSize = size(spekLabelBase,2); % number of x positions, this is also length(spekLabelBase)

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

numRand = length(testResults{iniCell}{1}.mvrgRand); % number of randomization done & shown in mvrg result

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

trajLabel = {'Immobile','Confined','Free','Directed','All'};
pSigTtest = NaN(groupsSize,numSmLabel,length(trajLabel));

%% Output
belowThreshFlag3D = false(groupsSize, numSmLabel, length(trajLabel));
fractionFlagMat3D = nan(groupsSize, numSmLabel, length(trajLabel));
avgPval = nan(groupsSize, numSmLabel, length(trajLabel));

%% Looping through 5 groups of trajectories in order: immobile, confined, free, directed, all together
coefMat4D = [];
for iRand = 1:numRand
    for iTraj = 1:5
        
        % %%%%% Collect MVRG coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if iTraj < 5 % section for different diffusion types
            
            coefMat3D = [];
            for iCell = 1 : numCell
                if ~isfield(testResults{iCell}{2}.mvrgs{iTraj},'nanFlag') && isstruct(testResults{iCell}{2}.mvrgs{iTraj})
                    coefMat3D = cat(3,coefMat3D,testResults{iCell}{2}.mvrgRand{iTraj}{iRand}.coef);
                else
                    missingData = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                    coefMat3D = cat(3,coefMat3D,missingData);
                end
            end % (for iCell = 1 : numCell)
 
        else % for all trajectories together (regardless of diffusion type)
              
            coefMat3D = [];
            for iCell = 1 : numCell
                if isstruct(testResults{iCell}{1}.mvrgs) && ~isfield(testResults{iCell}{1}.mvrgs,'nanFlag')
                    coefMat3D = cat(3,coefMat3D,testResults{iCell}{1}.mvrgRand{iRand}.coef);
                else
                    missingData = NaN(groupsSize,numSmLabel);%takes care of missing data in case mvrg is not run
                    coefMat3D = cat(3,coefMat3D,missingData);
                end
            end % (for iCell = 1 : numCell)
 
        end % (if iTraj < 5)
        
        
        % %%%%% Loop through each single-molecule property, perform outlier detection & ttest %%%%%%%%%%%%%%%%%%%%%%%%%%
        for iSMProp = 1 : size(coefMat3D,2)
            
            % Collect relevant coefficients
            tmp = squeeze(coefMat3D(1:end,iSMProp,:))';
 
            % Outlier analysis
            tmpFlag = false(size(tmp,1),groupsSize); % initialize tmpFlag
            for icol = 1:groupsSize
                
                % remove NaN before outlier detection, because sometimes isoutlier can be buggy when 10% of array > number of non-NaN data
                tmpCurrCol = tmp(:,icol); 
                indxNoNan = find(~isnan(tmpCurrCol));
                
                tmpFlag(indxNoNan,icol) = isoutlier (tmpCurrCol(indxNoNan),'gesd'); % detecting outliers
                %tmpFlag(:,icol) = isoutlier (tmp(:,icol),'gesd');% detecting outliers
                
            end % (for icol = 1:groupsSize)
            
            outlierIndxtmp = (tmpFlag == 1); % get linear indx of outlier,to makeoutliers NaN

            tmp(outlierIndxtmp) = NaN; % NaN the outliers in tmp 
            
            [~,pValue] = ttest(tmp, 0, 'Alpha', alphaTtest); % pvalue for only the inliers; (0 = null hypothesized population mean)
            
            pSigTtest(:,iSMProp,iTraj) = pValue';
            
            
        end % (for iSMProp = 1 : size(coefMat3D,2))
        
    end % (for iTraj = 1:5)
    
    coefMat4D = cat(4, coefMat4D, pSigTtest); % p-value from ttest of each randomization cases
    
end % (for iRand = 1:numRand)

sigMatFlag4D = coefMat4D <= alphaTtest; % indicates whether each randomization case contain significant coefficients by ttest with alpha-value = alphaTtest

fractionFlagMat3D = sum(sigMatFlag4D, 4)/numRand; % indicates fraction of randomization (0 to 1) contains significant coefficients by ttest with alpha-value = alphaTtest

belowThreshFlag3D = fractionFlagMat3D < alphaRandTest; % indicates whether each coefficients is significant by randomization test with alpha-value = alphaRandTest

pValSigOnly4D = coefMat4D.*sigMatFlag4D; % p-value from ttest of randomization cases that are deemed significant with alpha-value = alphaTtest
pValSigOnly4D(pValSigOnly4D == 0) = NaN; % Set Zeros To ‘NaN’

avgPval = mean(pValSigOnly4D, 4, 'omitnan'); % average of p-value from ttest of randomization cases that are deemed significant with alpha-value = alphaTtest
end