function plotBiphasicLines(TRI, motionType, propVect, threshTab, cellNo, figStr, scalingFactor)
%PLOTBIPHASICLINES: plot the fit lines for an individual cell from threshold and scatter data
%
%SYNOPSIS: plotBiphasicLines(TRI, motionType, propVect, threshTab, cellNo, figStr, scalingFactor)
%
%EXAMPLE: plotBiphasicLines(TRI, 2, [3 11], threshTab, 5, 'free')
%
%INPUT:    TRI        : Data of individual SM-speckle property to be
%                       plotted, can be input in either of the two following formats:
%
%                       (n x 1 cell array) MVRG result for individual
%                       cells. Must contain the field .aggMat for different
%                       motion types (immobile, confined, free, directed).
%                       Example: TRI = strInl.testResultAndSplit{1, 1};
%                       See runSmactinAggMovExplain.m for details.
%                                   ----------- OR -----------
%                       (nx3 matrix) 1st col = x-Data, 2nd col = y-Data,
%                       3rd col row 1 is threshold for this manual dataset.
%                   
%
%           motionType: (interger) motion type to pull data from (0 =
%                       immobile, 1 = confined, 2 = free, 3 = directed, 4 =
%                       all).
%
%           propVect  : (1 x 2 vector) columns of .aggMat that are to be plotted
%                     propVect(1) = x-axis data (e.g. col 3 = speckle density).
%                     propVect(2) = y-axis data (e.g. col 11= sm diffusion
%                                   coefficients, when there are 8 speckle properties).
%                       See smactinAggMov.m for details
%
%           threshTab : (table) threshold from biphasic-scanning, must
%                      include columns named:
%                           .thresholds
%                           .coefsOfMin
%                       See scan4smSpeckleThreshold.m for details.
%
%           cellNo    : (interger) cell from TRI (nx1 cell array) that is
%                       to be plotted. This is also the row of table
%                       threshTab from which the thresholds and fit 
%                       coefficients will be collected.
%
%           figStr    : (string) title of output figure
%
%           scalingFactor: (1x2 matrix) optional. Default [] = no scaling
%                       scalingFactor(1) = scaling for x-axis
%                       scalingFactor(2) = scaling for y-axis
%
%OUTPUT: figure plotting the fit lines for a particular cell on top of scatter data
%        2-line: red; 1-line: black.
%
% Tra Ngo, Oct 2021
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

%% Input
if ~exist('scalingFactor','var')
    scalingFactor = [];
end

%% Get x & y data and their range:
if iscell(TRI)
    if motionType <= 3
        xDat = TRI{cellNo,1}{2, 1}.aggMat{1, motionType+1}(:,propVect(1));
        yDat = TRI{cellNo,1}{2, 1}.aggMat{1, motionType+1}(:,propVect(2));
    elseif motionType == 4 % get data for all motions
        xDat = TRI{cellNo,1}{1, 1}.aggMat(:,propVect(1));
        yDat = TRI{cellNo,1}{1, 1}.aggMat(:,propVect(2));
    else
        error('Please input correct motionType (integer from 0 to 4)')
    end
    
    thres = threshTab.thresholds(cellNo);
    
else
    xDat = TRI(:,1);
    yDat = TRI(:,2);
    thres = TRI(1,3);
    
end % (if iscell(TRI))

x = min(xDat):(max(xDat)-min(xDat))/100:max(xDat);
x2lineBelow = min(xDat):(thres-min(xDat))/100:thres;
x2lineAbove = thres:(max(xDat)-thres)/100:max(xDat);

if ~isempty(scalingFactor)
    xDat = xDat*scalingFactor(1);
    yDat = yDat*scalingFactor(2);
    thres = thres*scalingFactor(1); %thres belongs to the x-axis
end

% Plot scatter and plot lines:
figure; scatter(xDat,yDat); hold on

% If there is no valid threshold, plot only scatter, no lines.
if isempty(threshTab.coefsOfMin{cellNo}) || any(any(isnan(threshTab.coefsOfMin{cellNo})))
    title(['Cell ' num2str(cellNo) ' ' figStr ' no threshold'])
    return
end

y1line = threshTab.coefsOfMin{cellNo}(1,3) + threshTab.coefsOfMin{cellNo}(2,3)*x;
yBelow = threshTab.coefsOfMin{cellNo}(1,1) + threshTab.coefsOfMin{cellNo}(2,1)*x2lineBelow;
yAbove =  threshTab.coefsOfMin{cellNo}(1,2) +  threshTab.coefsOfMin{cellNo}(2,2)*x2lineAbove;

if ~isempty(scalingFactor)
    y1line = y1line*scalingFactor(2);
    yBelow = yBelow*scalingFactor(2);
    yAbove = yAbove*scalingFactor(2);
    
    x = x*scalingFactor(1);
    x2lineBelow = x2lineBelow*scalingFactor(1);
    x2lineAbove = x2lineAbove*scalingFactor(1);
end

plot(x2lineBelow,yBelow,'r') % plot line below
plot(x2lineAbove,yAbove,'r') % plot line above
plot(x,y1line,'k') % plot 1 line

xline(thres,'--r') % plot biphasic threshold

title(['Cell ' num2str(cellNo) ' ' figStr])

end