function [outlierLabelVect] = labelOutliersFrMat(dataMat, neighborNum, k, exclCol)
%LABELOUTLIERSFRMAT  reads n-rows of observations with m-cols of propertiese and label which observation (row) is an outlier
%
%SYNOPSIS:  [outlierLabelVect] = labelOutliersFrMat(dataMat, neighborNum)
%
%INPUT:     dataMat: (matrix) data that need outliers detected
%                   Each row: a data point (an observation)
%                   Each colum: a feature (a dimesion of the data point)
%
%           neighborNum: (integer) number of nearest neighbor to calculate
%                   distance of a particular observation from.
%                   TN20210615: remove option to leave this as default = 10.
%
%           k : (optional, integer) Roughly, for a certain value of k, 
%               observations that are k*sigma away from the mean will be
%               considered outliers. (default: [], so that value is
%               determined by detectOutliers.m's default input).
%
%           exclCol: (optional, integer) column in dataMat where identical
%                    values will be searched for and removed before detection 
%                    of outliers. (default: no identical datapoint in any
%                    colume is searched for & removed) 
%
%OUTPUT:    outlierLabelVect: logic vector denoting at each index whether
%                       the dataVect's entry is an outlier (0) or not (1)
%
%REMARKS: for an n-by-m dataMat, functions operate in the following steps
%STEP 1: calculate distance (L2 norm) to 10 nearest neighbor & average the mean
%STEP 2: use detectOutlier.m on this 1D data
%STEP 3: label identified outliers (a new logic vector: 1 = in-lier & 0 = outliers)
%
%Tra Ngo, Oct 2020
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

%% Input control
if ~exist('neighborNum','var') || isempty(neighborNum)
    error('Please input the number of neighbor for averaging pairwise distance.');
elseif neighborNum > size(dataMat,1)
    warning('Number of neighbor is larger than number of available datapoint.');
end
if ~exist('k','var')
    k = [];
end

[nRow, ~] = size(dataMat);

datExcludeFlag = ~exist('exclCol', 'var') || isempty(exclCol) || exclCol == 0; % TN20210629 Default to not exclude anything in dataMat

%% Initiate variables, constants, and output
outlierLabelVect = true(nRow,1); % note that logic has much lower storage requirement than ones()
figFlag = false; % for code development and debugging purpose


%% STEP1: calculate distance to nearest neighbors
pairwiseDist = pdist(dataMat); % compare speed against createDistanceMatrix.m


if ~datExcludeFlag % remove duplicate datapoint
    % TN20210608: pdist of diffusion coef alone, find index of 0 (indicates distance
    % between 2 datapoints is 0, meaning they are identical) => replace with NaN.
    exclDist = pdist(dataMat(:,exclCol));
    duplIndx = (exclDist == 0); % logical index of pairwise distance = 0.
    pairwiseDist(duplIndx) = NaN;
end

fullDistMat = squareform(pairwiseDist);

%checking distance between data#1 and all data-point & find k-smallest number of pDistRes
% (neighborNum+1) col because distance to self (=0) is included in fullDistMat 
distVectPerDataMat = mink(fullDistMat, (neighborNum + 1), 2); % n-by-(neighborNum+1) matrix

% take the average for each data's distance to nearest neighbors
meanDistVect = sum(distVectPerDataMat,2)/neighborNum; % n-by-1 vector

%% STEP 2: use detectOutlier.m on this 1D data
%[outlierIdx,~] = detectOutliers(meanDistVect, k);
%outlierIdx = isoutlier(meanDistVect,'gesd'); % 1 = outlier, 0 = inlier; 20210622 TN changed to GESD outlier detection
outlierIdx = isoutlier(meanDistVect,'quartiles','ThresholdFactor',k); % 1 = outlier, 0 = inlier; 20210624 TN changed to Quartiles outlier detection

%% STEP 3: label identified outliers (a new logic vector: 1 = in-lier & 0 = outliers)
outlierLabelVect(outlierIdx) = false;

%% Visualize mean distance and outliers
if figFlag % for code development and debugging purpose
    figure
    histogram(meanDistVect); hold on, xlabel('mean distance to nearest neighbors')
    try scatter(meanDistVect(outlierIdx),3*ones(1,sum(outlierIdx)));
    catch, scatter(meanDistVect(outlierIdx),3*ones(1,length(outlierIdx))); end
    legend('histogram','Outliers'), title('Histogram of mean distance with outliers denoted')
    

    figure
    anscombMeanDistVect = sqrt(meanDistVect + 3/8);
    histogram(anscombMeanDistVect); hold on, xlabel('Anscombe-transformed mean distance to nearest neighbors')
    try scatter(anscombMeanDistVect(outlierIdx),3*ones(1,sum(outlierIdx)));
    catch, scatter(anscombMeanDistVect(outlierIdx),3*ones(1,length(outlierIdx))); end
    legend('histogram','Outliers'), title('Anscombe-transformed 1D distance')
end

end