function [scoreTot, frameIndexTot, meanScoreByFrmMat] = speckleFrameAnalyse(actinPropPerTimeInt, outPath, fieldName, plotType, splitNum, figName)
%SPECKLEFRAMEANALYSE  plots speckle properties for each frame of all movie (i.e. kinetic scores, speckle speed, etc) in an actinPropPerTimeInt structure
%
%SYNOPSIS   [scoreTot, frameIndexTot, meanScoreByFrmMat] = speckleFrameAnalyse(actinPropPerTimeInt,outPath,fieldName,plotType,splitNum, figName)
%
%EXAMPLE: speckleFrameAnalyse(actinPropPerTimeInt,'/project/biophysics/jaqaman_lab/vegf_tsp1/...','kinScore','box',1);
%
%INPUT: actinPropPerTimeInt: (n x 1 cell array of structure) speckles
%                               information from FSM movies. See
%                               actinPropPerTimeInterval.m for details.
%
%       outPath            : (string) folder where generated plots are to be saved.
%
%       fieldName          : (string) name of the field that contain the
%                            information to be plotted, can be one of the following:
%                                  .kinScorePos
%                                  .kinScore
%                                  .speckleList
%                                  .speckleInitPos
%                                  .speckleMidPos
%                                  .speckleVelocity
%                                  .speckleSpeed
%                                  .speckleDensity
%                                  .speckleMvmtCohere
%                                  .indxMask
%                                  .speckleMaskDensity
%
%       plotType           : (string) indicate the kind of plot to visualize 
%                            data in, can be one of the following strings:
%                                          'box'
%                                          'stair'
%                                          'notBoxPlot'
%
%       splitNum           : (integer) number of frame apart to be plotted.
%
%       figName            : (string) added to the end of figure name.
%
%OUTPUT: generated plots will be saved in the outPath indicated from the
%           input.
%
%        scoreTot          : (array) property of interest that is plotted.
% 
%        frameIndexTot     : (vector) the frame which corresponding
%                            property is pulled from.
% 
%        meanScoreByFrmMat : (matrix) means of scoreTot.
%
%REMARK     Note that there is an assumption that the FSM movies are 19 frames.
%
%Tra H. Ngo, February 2019
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

%% definition
frameStart = 1;
try
    frameEnd = size(actinPropPerTimeInt{1},1);
catch
    
    frameEnd = 19; % this function only accomodates movies of 19 frames right now.
end
binSize = 0.0002; % this number is just a place holder and will be empirically defined for each case.
movNum = size(actinPropPerTimeInt,1);
switch fieldName
    case 'kinScore'
        figTitle = 'Kinetic Score';
        frameStart = 2; % we don't have any kinetics scores for first frame (1) and last frame (19) (which are 0) as we cannot calculate any derivatives for them
        binSize = 0.0002; % for making figure purpose, this value is empirical since our kinetic score usually doesn't exceed abs(0.008)
    case 'speckleSpeed'
        figTitle = 'Speckle Speed';
        binSize = 0.2; % for making figure purpose, this value is empirical since our speckle speeds usually differ by 0.2
    case 'speckleDensity'
        figTitle = 'Speckle Density';
        binSize = 0.5; % for making figure purpose, this value is empirical since our speckle speeds are interger
    case 'speckleVelocity'
        figTitle = 'Speckle Velocity';
    case 'speckleMvmtCohere'
        figTitle = 'Speckle Movement';
        binSize = 0.2;
    case 'speckleMaskDensity'
        figTitle = 'Speckle Mask Density';
        binSize = 0.001;
    case 'ILmax'
        figTitle = 'Speckle Intensity';
        binSize = 0.001;
    otherwise
        figTitle = 'Default fig name';
        binSize = 0.001;
        warning('Plots may not display name of appropriate requested field or appropriate binSize.')
end

%% output initialization
scoreTot = []; %initiate matrix of total score
frameIndexTot = []; %initiate matrix of total kinetic score
meanScoreByFrmMat = [];

%% Plot and save different kinds of plot
switch plotType
    case 'stair'
        %% Stair plot for each frame of all movie
        maxNum2 = max(abs(actinPropPerTimeInt{1}(2).(fieldName)));
        maxNum3 = max(abs(actinPropPerTimeInt{1}(3).(fieldName)));
        maxNum4 = max(abs(actinPropPerTimeInt{1}(4).(fieldName)));
        % to estimate the range of the values that are being ploted, we take the
        % maximum of the biggest value from randomly-selected frame #2,3,4 of the first movie
        figLim = round(max([maxNum2,maxNum3,maxNum4]),4);
        binEdges = (-(figLim + binSize):binSize:(figLim + binSize));
        
        % binNum = length(binEdges)-2;
        
        frameCol = linspace(0.3,1,21);
        f = figure;
        hold on
        
        for iFrame = frameStart:splitNum:frameEnd
            countPos = zeros;
            for iMov = 1:movNum
                countPos = countPos + histcounts(actinPropPerTimeInt{iMov}(iFrame).(fieldName),binEdges);
            end
            stairPosNorm = countPos/sum(countPos);
            
            stairs(binEdges(1:end-1),stairPosNorm, 'DisplayName', ...
                ['Frame', num2str(iFrame), ' with ', num2str(sum(countPos)),' ', figTitle] , ...
                'color',[0.5 frameCol(iFrame) frameCol(iFrame)])
        end
        legend('show')
        title(['All Movies All Frames Normalized ', figTitle, ' ', figName]);
        
        figureName = [fieldName, '_Normalized_' figName '.fig'];
        
        
        
    case {'box', 'violin', 'notBoxPlot','line'}
        %% Box plot
        
        for iMov = 1:movNum
            if isempty(actinPropPerTimeInt{iMov})
                continue
            end
                
            switch fieldName
                case 'kinScore' % last frame does not have any kinScore value
                    frameEnd = length(actinPropPerTimeInt{iMov}) - 1;
                otherwise
                    frameEnd = length(actinPropPerTimeInt{iMov});
            end
            
            ScoreCell = vertcat(actinPropPerTimeInt{iMov}(frameStart:frameEnd).(fieldName));
            
            
            for iFrame = frameStart:frameEnd
                % create a vector indicating which frame a corresponding score is pulled from
                frameIndex = iFrame*ones(length(actinPropPerTimeInt{iMov}(iFrame).(fieldName)),1);
                frameIndexTot = vertcat(frameIndexTot, frameIndex);
                
                % get mean value of each frame, for each cell
                switch plotType
                    case {'notBoxPlot','line'}
                        meanScoreByFrm{iMov,1}(iFrame) = nanmean(actinPropPerTimeInt{iMov}(iFrame).(fieldName));
                    otherwise
                end
            end
            
            scoreTot = vertcat(scoreTot, ScoreCell);
            
        end % (for iMov = 1:movNum)
        
        
        % Plot
        f = figure;
        
        switch plotType
            case 'box'
                boxplot(scoreTot, frameIndexTot);
                
                figureName = [fieldName, '_BoxPlot_' figName '.fig'];
                
            case 'notBoxPlot'
                meanScoreByFrmMat = cell2mat(meanScoreByFrm); % row = cell, col = frame
                notBoxPlot(meanScoreByFrmMat, 1:size(meanScoreByFrmMat,2),[],[],0)
                
                figureName = [fieldName, '_notBoxPlot_' figName '.fig'];
                
            case 'line'
                meanScoreByFrmMat = cell2mat(meanScoreByFrm); % row = cell, col = frame
                plot(meanScoreByFrmMat', '-o')
                
                figureName = [fieldName, '_linePlot_' figName '.fig'];
                
            case 'violin'
                switch fieldName  % consistent speckleSpeed through the frames will show ksdensity error!
                    case 'speckleSpeed'
                        violinplot(scoreTot, frameIndexTot);
                    otherwise
                        violinplot(scoreTot, frameIndexTot, 'ShowNotches',true);
                end
                
                figureName = [fieldName, '_ViolinPlot_' figName '.fig'];
                
            otherwise
        end % (switch plotType)
        
        xlabel('Frame #')
        title(['All Movies All Frames: ', figTitle, ' ', figName]);
        
        
    otherwise
        warning('Function does not recognize this type of plot.')
end % (switch plotType)

%% Save box plot in output folder
cd(outPath) %change directory to path to look for folder
saveas(f, figureName);

end