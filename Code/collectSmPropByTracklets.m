function [propPerSmPerCellOutput, propStatsPerCellOutput, trackClassSizePerCellOutput, inputParam] = ...
    collectSmPropByTracklets(combnPropPerTimeInt, prop2extract, indicesOfInt, figFlag)
%COLLECTSMPROPBYTRACKLETS: collect specified single-molecule property from combinePropPerTimeInt
%
%SYNOPSIS:  [propPerSmPerCellOutput,propStatsPerCellOutput,trackClassSizePerCellOutput, inputParam] = ...
%                   collectSmPropByTracklets(combnPropPerTimeInt,prop2extract,indicesOfInt, figFlag)
%
%INPUT:    combnPropPerTimeInt : Cell array with length = number of movies.  
%                                 See smactinPropPerTimeInt.m for details.
%
%          prop2extract        : (char.) name of field of SM property that
%                           user wants to collect from combinePropPerTimeInt.
%
%          indicesOfInt        : (vector of integer) number denoting which 
%                           cells from cell array of combinePropPerTimeInt
%                           to collect information from.
%
%          figFlag             : (logical) flag denoting whether to output figures.
%                           1 = (Default) output and save figures in a new subfolder
%                           0 = not output figures.
%
%OUTPUT:   (with n being number of movies input in combnPropPerTimeInt)
%
%          propPerSmPerCellOutput: (n x 5 x 3 cells):
%
%                        1st cell: (n x 5 cell) n row for each cell, 5 columns: total 
%                                  and then the 4 diffusion classes; each entry 
%                                  contains the SM property of interest. SM 
%                                  regardless of number of matched speckle.
%                        2nd cell: (n x 5 cell) similar to 1st cell but for
%                                  SM matched to 0 speckle.
%                        3rd cell: (n x 5 cell) similar to 1st cell but
%                                  for SM matched to >=1 speckles.
%
%          propStatsPerCell    : (n x 5 x 3 cells)
%                        1st cell: SM info regardless of the number of matched speckle,
%                        2nd cell is SM info for SM matched to 0 speckle, 3rd cell is 
%                        SM info for SM matched to >=1 speckles. Each nx5 cell contains 
%                        n row for each cell, 5 columns: total and then the 4 diffusion
%                        classes; each entry contains 1x3 matrix of [mean, median, std].
%
%          trackClassSizePerCell: (1 x 3 cell)
%                        1st cell (nx7x3 matrix) contains raw number of SM track (the
%                                 columns are: total & 4 diffusion classes & non-classified 
%                                 tracks & total classified tracks) for n movies.
%                        2nd cell (nx4x3 matrix)contains fraction over total classified  
%                                 tracks (only the 4 diffusion classes)); n rows with 
%                                 each row is one movie. 
%                        3rd cell table summing up all tracks of each group in 1st cell.
%
%                        1st page is info regardless of number of matched speckle,
%                        2nd page is SM info for SM matched to 0 speckle,
%                        3rd page is SM info for SM matched >= 1 speckles.
%                        Previously 2nd cell contains fraction over total (only the 4 
%                        diffusion classes), while 3rd cell contains fraction over 
%                        classified tracks. (TN 20210901)
%
%          inputParam          : = {prop2extract,indicesOfInt}
%
%Tra Ngo, Aug. 2021 (modified from Aparajita's script)
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

%% Input control & variable initiation

numCell = length(combnPropPerTimeInt);
numIndx = length(indicesOfInt);

if numIndx > numCell
    error('Number of cells being asked to collect can only equal or smaller than the number of available cells.')
end

% Check if input fieldname is of type 'char'. This does not check if this
% field exist within each cell's structure.
if ~exist('prop2extract','var') || ~ischar(prop2extract) || isempty(prop2extract)
    error('Property to extract must be input as character variable type.')
end

if ~exist('figFlag','var') || ~islogical(figFlag)
    figFlag = true;
end

%% Output initiation & variable initiation
propPerSmPerCell = cell(numIndx, 5);
propPerSmNoSpkl = cell(numIndx, 5); % SM property collected here belongs to SM that do not get matched to any speckle
propPerSmYesSpkl = cell(numIndx, 5); % SM property collected here belongs to SM that do not get matched to multiple (1 or more) speckles
propPerSmPerCellOutput = cat(3,propPerSmPerCell,propPerSmNoSpkl,propPerSmYesSpkl);

propStatsPerCell = cell(numIndx, 5);
propStatsNoSpkl = cell(numIndx, 5);
propStatsYesSpkl = cell(numIndx, 5);
propStatsPerCellOutput = cat(3, propStatsPerCell, propStatsNoSpkl, propStatsYesSpkl);

trackClassNumPerCell = nan(numIndx, 7); % total and the 4 diffusion classes & non-classified & total classified tracks
trackClassNumNoSpkl = nan(numIndx, 7); % total and the 4 diffusion classes & non-classified & total classified tracks
trackClassNumYesSpkl = nan(numIndx, 7); % total and the 4 diffusion classes & non-classified & total classified tracks
trackClassFracPerCell = nan(numIndx, 4); % only the 4 diffusion classes
trackClassFracNoSpkl = nan(numIndx, 4); % only the 4 diffusion classes, percentage matched to 0 speckle
trackClassFracYesSpkl = nan(numIndx, 4); % only the 4 diffusion classes, percentage matched to >= 1 speckles

inputParam = {prop2extract, indicesOfInt};

if numCell == 0 || numIndx == 0
    warning('User input no cell to collect. Function returns to main.')
    return
end


%% Get SM information
for i = 1:numIndx % loop through the cells that user wants to collect info from
    iMov = indicesOfInt(i);
    if isempty(combnPropPerTimeInt{iMov, 1})
        continue
    end
    
    numFrame = length(combnPropPerTimeInt{iMov, 1}{1, 2}); % number of FSM frame in this cell
    
    
    classVectAll = []; propVectAll = []; numMatchVctAll = [];
    for iFr = 1:numFrame % loop throught the FSM frames of each cell
        
        % getting trackclass info
        classVect = combnPropPerTimeInt{iMov, 1}{1, 2}(iFr).trackClass;
        classVectAll = vertcat(classVectAll, classVect); % vertcat of classVect from all FSM frame
        % getting property of interest
        propVect = combnPropPerTimeInt{iMov, 1}{1, 2}(iFr).(prop2extract); % all diffusion type
        propVectAll = vertcat(propVectAll, propVect);
        
        % getting list of matched speckle
        spkList = combnPropPerTimeInt{iMov, 1}{1, 2}(iFr).speckleList; % all diffusion type
        if isempty(spkList) % account for when first FSM interval does not have corresponding SM information (LatA 2018 movies)
            continue
        end
        numMatch = cellfun(@length, spkList, 'UniformOutput',false);
        numMatchVct = cell2mat(numMatch);
        numMatchVctAll = vertcat(numMatchVctAll, numMatchVct);
        
    end % (for iFr = 1:numFrame)
    
    propPerSmPerCell{i,1} = propVectAll; % all diffusion type, all frames
    
    % getting statistic of SM's property
    meanTmp = nanmean(propPerSmPerCell{i,1}); 
    medianTmp = nanmedian(propPerSmPerCell{i,1}); 
    stdTmp = nanstd(propPerSmPerCell{i,1}); 
    propStatsPerCell{i,1} = [meanTmp medianTmp stdTmp]; 
    
    % getting total number of SM
    trackClassNumPerCell(i,1) = length(classVectAll);
    
    % getting number of unclassified SM tracks
    indxNan = isnan(classVectAll);
    trackClassNumPerCell(i,6) = sum(indxNan);
    % getting number of classified SM tracks
    trackClassNumPerCell(i,7) = trackClassNumPerCell(i,1) - trackClassNumPerCell(i,6);
    
    % Get index SM with corresponding number of matched speckles
    indxMtchNone = numMatchVctAll == 0;
    indxMtchMultpl = numMatchVctAll >= 1;
    
    % Get property of SM that get matched to 0 or >= 1 speckle
    propPerSmNoSpkl{i,1} = propPerSmPerCell{i,1}(indxMtchNone,:);
    propPerSmYesSpkl{i,1} = propPerSmPerCell{i,1}(indxMtchMultpl,:);
    
    % getting statistic of SM's property when there are 0 or >=1 matched speckle(s)
    meanTmp = nanmean(propPerSmNoSpkl{i,1}); 
    medianTmp = nanmedian(propPerSmNoSpkl{i,1}); 
    stdTmp = nanstd(propPerSmNoSpkl{i,1}); 
    propStatsNoSpkl{i,1} = [meanTmp medianTmp stdTmp];
    
    meanTmp = nanmean(propPerSmYesSpkl{i,1}); 
    medianTmp = nanmedian(propPerSmYesSpkl{i,1}); 
    stdTmp = nanstd(propPerSmYesSpkl{i,1}); 
    propStatsYesSpkl{i,1} = [meanTmp medianTmp stdTmp];
    
    % getting number of SM tracks with and without speckle
    trackClassNumNoSpkl(i,1) = sum(indxMtchNone);
    trackClassNumYesSpkl(i,1) = sum(indxMtchMultpl);
    
   for iType = [0 1 2 3]
        
        % Store information of SM from motion type 0 immobile in col 2,
        % type 1 confined in col 3, type 2 free in col 4, type 3 directed in col 5.
        indexType = classVectAll == iType;
        
        propPerSmPerCell{i,iType+2} = propPerSmPerCell{i,1}(indexType,:);
        
        % Store statistic of SM's property
        meanTmp = nanmean(propPerSmPerCell{i,iType+2});
        medianTmp = nanmedian(propPerSmPerCell{i,iType+2});
        stdTmp = nanstd(propPerSmPerCell{i,iType+2});
        propStatsPerCell{i,iType+2} = [meanTmp medianTmp stdTmp];
        
        % Store number of SM per trajectory
        numTrackTmp = sum(indexType); % calling "sum" is faster than "length" here.
        trackClassNumPerCell(i,iType+2) = numTrackTmp;
        trackClassFracPerCell(i,iType+1) = trackClassNumPerCell(i,iType+2)/...
            trackClassNumPerCell(i,7); % divide by total number of CLASSIFIED tracks
        
       % Get property of SM that get matched to 0 or >= 1 speckle
        propPerSmNoSpkl{i,iType+2} = propPerSmPerCell{i,iType+2}(indxMtchNone(indexType),:);
        propPerSmYesSpkl{i,iType+2} = propPerSmPerCell{i,iType+2}(indxMtchMultpl(indexType),:);
        
        % Store statistic of SM's property when there are 0 or >=1 matched speckle(s)
        meanTmp = nanmean(propPerSmNoSpkl{i,iType+2});
        medianTmp = nanmedian(propPerSmNoSpkl{i,iType+2});
        stdTmp = nanstd(propPerSmNoSpkl{i,iType+2});
        propStatsNoSpkl{i,iType+2} = [meanTmp medianTmp stdTmp];
        
        meanTmp = nanmean(propPerSmYesSpkl{i,iType+2});
        medianTmp = nanmedian(propPerSmYesSpkl{i,iType+2});
        stdTmp = nanstd(propPerSmYesSpkl{i,iType+2});
        propStatsYesSpkl{i,iType+2} = [meanTmp medianTmp stdTmp];
        
        % Store number of SM per trajectory that get matched to different number of speckles
        trackClassNumNoSpkl(i,iType+2) = sum(indxMtchNone(indexType)); % matched to 0 speckle
        trackClassNumYesSpkl(i,iType+2) = sum(indxMtchMultpl(indexType)); % matched to >= 1 speckles
        
    end % (for iType = [0 1 2 3])
    
    trackClassNumNoSpkl(i,7) = sum(trackClassNumNoSpkl(i,2:5)); % total classified
    trackClassNumNoSpkl(i,6) = trackClassNumNoSpkl(i,1) - trackClassNumNoSpkl(i,7); % total unclassified
    trackClassNumYesSpkl(i,7) = sum(trackClassNumYesSpkl(i,2:5)); % total classified
    trackClassNumYesSpkl(i,6) = trackClassNumYesSpkl(i,1) - trackClassNumYesSpkl(i,7); % total unclassified
    
    for iType = [0 1 2 3]
        trackClassFracNoSpkl(i,iType+1) = trackClassNumNoSpkl(i,iType+2)/trackClassNumNoSpkl(i,7); % divide by total number of CLASSIFIED tracks matched to 0 speckle
        trackClassFracYesSpkl(i,iType+1) = trackClassNumYesSpkl(i,iType+2)/trackClassNumYesSpkl(i,7);% divide by total number of CLASSIFIED tracks matched to >= 1 speckle
    end % (for iType = [0 1 2 3])
    
end % (for i = 1:numIndx)

totNumTrack = [sum(trackClassNumPerCell,1); sum(trackClassNumNoSpkl,1); sum(trackClassNumYesSpkl,1)];
rowNames = {'All tracks','0 speckle','> 0 speckles'};
colNames = {'All','Immobile','Confined','Free','Directed','Unclassified', 'Classified'};
sTable = array2table(totNumTrack,'RowNames',rowNames,'VariableNames',colNames);

if figFlag
    %% Plot not box plot comparing mean of immobile, confined, free, directed
    tmp = cell2mat(propStatsPerCell);
    h = figure;
    notBoxPlot(tmp(:,1:3:15),[],[],[],0) % tmp(:,[1:3:15]) contains 5 col, each for all & motion type; each row is the mean property value of each cell (movie)
    title(['Mean of ' prop2extract ' (each dot is a cell)']);
    h.Children.XTickLabel = {'all', 'immobile', 'confined', 'free', 'directed'};
    
    %% Plot not box plot comparing mean of immobile, confined, free, directed
    tmp = cell2mat(propStatsNoSpkl);
    h2 = figure;
    notBoxPlot(tmp(:,1:3:15),[],[],[],0) % tmp(:,[1:3:15]) contains 5 col, each for all & motion type; each row is the mean property value of each cell (movie)
    title(['Mean of ' prop2extract ' (each dot is a cell), no speckle']);
    h2.Children.XTickLabel = {'all', 'immobile', 'confined', 'free', 'directed'};
    
    %% Plot not box plot comparing mean of immobile, confined, free, directed
    tmp = cell2mat(propStatsYesSpkl);
    h3 = figure;
    notBoxPlot(tmp(:,1:3:15),[],[],[],0) % tmp(:,[1:3:15]) contains 5 col, each for all & motion type; each row is the mean property value of each cell (movie)
    title(['Mean of ' prop2extract ' (each dot is a cell), >0 speckles']);
    h3.Children.XTickLabel = {'all', 'immobile', 'confined', 'free', 'directed'};
    
    %% Plot violin plot containing aggregated properties of interest of all cells
    propEach = cell(5,1); lenDat = nan(5,1);
    for j = 1:5
        propEach{j} = cell2mat(propPerSmPerCell(:,j));
        lenDat(j) = length(propEach{j});
    end
    propTmp = padcat(propEach{1}, propEach{2}, propEach{3}, propEach{4}, propEach{5});
    f = figure;
    title([prop2extract ' of each motion type, all cells']);
    violinplot(propTmp);
    f.Children.XTickLabel = {['all (' num2str(lenDat(1)) ')'], ['immobile (' num2str(lenDat(2)) ')'], ...
        ['confined (' num2str(lenDat(3)) ')'], ['free (' num2str(lenDat(4)) ')'], ...
        ['directed (' num2str(lenDat(5)) ')']};
    
    
    %% Plot notBoxPlot to show fraction of motion types among only classified receptors
    f3 = figure;
    title('Fraction of motion types among classified receptors')
    %notBoxPlot(trackClassSizePerCell{1,3},[],[],[],0)
    notBoxPlot(trackClassFracPerCell,[],[],[],0)
    f3.Children.XTickLabel = {'immobile', 'confined', 'free', 'directed'};
    
    %% Plot notBoxPlot to show fraction of motion types among those matched to no speckle
    f4 = figure;
    title('Fraction of receptors matched to no speckle among each motion type')
    notBoxPlot(trackClassFracNoSpkl,[],[],[],0)
    f4.Children.XTickLabel = {'immobile', 'confined', 'free', 'directed'};
    
    %% Plot notBoxPlot to show fraction of motion types among those matched to >= 1 speckles
    f5 = figure;
    title('Fraction of receptors matched to one or more speckles among each motion type')
    notBoxPlot(trackClassFracYesSpkl,[],[],[],0)
    f5.Children.XTickLabel = {'immobile', 'confined', 'free', 'directed'};
    
    %% save all figure in current directory:
    tmp = clock; tmp = num2str(tmp(1:5)); tmp(tmp == ' ') = [];
    newFolder = [tmp prop2extract 'Figures'];
    mkdir(newFolder); cd(newFolder)
    saveas(h,'meanOfProp.fig')
    saveas(h2,'meanOfProp_noSpeckle.fig')
    saveas(h3,'meanOfProp_yesSpeckle.fig')
    saveas(f,'violinOfProp.fig')
    saveas(f3, 'receptorOfFraction_total.fig')
    saveas(f4, 'receptorOfFraction_noSpeckle.fig')
    saveas(f5, 'receptorOfFraction_yesSpeckle.fig')
    cd ..
    
end

% Final output
propPerSmPerCellOutput = cat(3,propPerSmPerCell,propPerSmNoSpkl,propPerSmYesSpkl);
propStatsPerCellOutput = cat(3, propStatsPerCell, propStatsNoSpkl, propStatsYesSpkl);
trackClassSizePerCellOutput = {cat(3, trackClassNumPerCell, trackClassNumNoSpkl, trackClassNumYesSpkl) cat(3,trackClassFracPerCell,trackClassFracNoSpkl,trackClassFracYesSpkl) sTable};

end