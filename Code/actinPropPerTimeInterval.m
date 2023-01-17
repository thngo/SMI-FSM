function [actinPropPerTimeInt] = actinPropPerTimeInterval(qfsmPackPaths, actinConditions)
%ACTINPROPPERTIMEINTERVAL : reads results from QFSM Package and calculates spatiotemporal properties of speckles per FSM frame
%
%SYNOPSIS: [actinPropPerTimeInt] = actinPropPerTimeInterval(qfsmPackPaths, actinConditions)
%
%INPUT:   qfsmPackPaths : (n x 1 cell) paths leading to QFSM Package results.
%
%       actinConditions : (structure) parameters for processing speckles;
%                         contains the following fields: 
%
%                        .rad_density: (n x 1 vector of integer) search radius for 
%                                     neighboring speckles for individual speckle 
%                                     density calculation. Recommended value: 7
%
%                        .maskLoc    : (n x 1 vector of integer) flag indicating
%                                     location of cell ROI mask.
%                            = 0 : there exists a folder with external
%                                  masks within path provided in
%                                  qfsmPackPaths (make sure folder name
%                                  follows a name in extlroi variable.)
%                            = 1 : there exists a folder with refined
%                                  masks within the QFSMPackage folder
%                                  provided in qfsmPackPaths
%                            = 2 : there exists a folder with refined
%                                  masks within the SegmentationPackage
%                                  folder provided in qfsmPackPaths
%                            If [] or other values, code will crash.
%                            (TN.AD. 20211119)
%
%                        .pos_std    : (n x 1 vector of integer) flag for
%                                     adding noise into speckle position data.
%                                    = 1 : gaussian noise will be added to speckle
%                                          positions. Noise follows Normal ~ (0,(15nm/90).^2)
%                                    = 0 or any value not 1: No noise added.
%
%                        .minLft     : (n x 1 vector of integer) minimum speckle 
%                                  lifetime to use in analysis. Recommended
%                                  value: 2 (to remove "ghost" speckles)
%
%OUTPUT: actinPropPerTimeInt: (n x 1 cell array of structures) FSM arm
%                             analysis results of SMI-FSM analysis for n
%                             movies; each structure contain speckle
%                             properties for m FSM intervals, containing
%                             the following fields:
%
%               .speckleList      : List of speckles that exist in the frame
%
%               .kinScore         : The kinetic score from the QFSM package
%
%               .kinScorePos      : The positions at which we have kinetic scores
%                                   (units: pixels)
%
%               .speckleInitPos   : The speckle position for a particular FSM
%                                   frame (prior frame of the interval) (units: pixels)
%
%               .speckleMidPos    : The speckle midpoint position between 2
%                                   FSM frames of an FSM interval (units: pixels)
%
%               .speckleVelocity  : The speckle velocity (displacement) from 
%                                   a frame to the next (units: pixels/frame)
%
%               .speckleSpeed     : The speckle speed (displacement magnitude) 
%                                   from one FSM frame to the next (units: pixels/frame)
%
%               .speckleMvmtCohere: The average comovement between a speckle 
%                                   and its neighboring speckles.
%
%               .maskDensity      : The number of speckles within a cell ROI 
%                                   mask divided by the number of pixels
%                                   comprising the mask.
%
%               .indxMask         : Indices of the mask obtained from QFSM package
%
%               .ILmax            : intensity of a specific speckle, as 
%                                   reported from QFSM package.
%
%               .IBkg             : background intensity around a specific
%                                   speckle, as reported from QFSM package.
%
%               .ILmaxNeighb      : average intensity of speckles in the
%                                   neighborhood of a specific speckle.
%
%               .IBkgNeighb       : average background intensity of speckles
%                                   in the neighborhood of a specific speckle.
%
% Deryl Tschoerner, February 2018
% modified: Tra H Ngo, February 2019 (remove ghost speckle)
% modified: Tra H Ngo, Aug. 2019 (add speckleLifetime properties)
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

%% definitions
numMovs = length(qfsmPackPaths);
fs = filesep;
maps = cell(numMovs,1);
actinPropPerTimeInt = cell(numMovs,1);
fsmFrames = nan(numMovs,1);

% %naming possible combinations of external masking names which were used in
% %the past. Fourth option may be applied most often as of June 2018
extlroi = ["external_roi","external_mask1","External_mask","External_ROI","external_ROI", ...
    "External_Segmentation"];
% as of February 2019, we have decided that refined_masks from Segmentation
% packages are the same as those from External_ROI, making option 1
% and 3 similar

for i = 1 : numMovs
    
    %field names for the variables found in the kineticAnalysis folder
    maps{i,1} = dir(strcat([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'kineticAnalysis' fs '*.mat']));
    fsmFrames(i,1) = length(maps{i,1});
    
    %skip an element in the cell if path does not exist.
    if isempty(qfsmPackPaths{i,1})
        continue
    end
    
    % Load speckle tracks (MPM)
    load([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'speckleTracks' fs '_tracks.mat'], 'MPM');
    
    
    % remove ghost speckles and calculate speckle lifetime properties
    [MPM,lifetimeMPM] = retainMPMminLifetime(MPM, actinConditions.minLft(i));
    
    %KJ 20210707: indicate complete vs. incomplete speckle lifetime
    flagCompleteLft = ones(size(lifetimeMPM));
    indxFirstFrame = find(lifetimeMPM(:,1) > 0);
    for iRow = indxFirstFrame'
        flagCompleteLft(iRow,1:2*lifetimeMPM(iRow,1)) = 0;
    end
    indxLastFrame = find(lifetimeMPM(:,end) > 0);
    for iRow = indxLastFrame'
        flagCompleteLft(iRow,end-2*lifetimeMPM(iRow,end)+1:end) = 0;
    end
    

    %TN2019: replace zeros with NaN
    MPM(MPM == 0) = NaN; MPM(:,end + 1 : end + 2) = NaN;
    
    if actinConditions.pos_std(i) == 1, noise = normrnd(0,15/90,size(MPM)); MPM = MPM + noise; end
    
    %store field names for readout of values from the structure fields found
    %in 'maps' variable
    temp = repmat(struct('polyMap',[],'kinScore',[],'depolyMap',[], ...
        'kinMap2C',[]),fsmFrames(i,1),1);
    
    %field names for vector of structures holding speckle properties
    actinPropPerTimeInt{i,1} = repmat(struct('kinScorePos',[],'kinScore',[], ...
        'speckleList',[],'speckleInitPos',[],'speckleMidPos',[], ...
        'speckleVelocity',[],'speckleSpeed',[],'speckleDensity',[] , ...
        'speckleMvmtCohere',[],'speckleLifetime',[],'flagCompleteLft',[], ...
        'ILmax', [], 'IBkg', [], 'ILmaxNeighb', [], 'IBkgNeighb', []),...
        fsmFrames(i,1),1);
    
    
    if ~isempty(qfsmPackPaths{i,1}) && actinConditions.maskLoc(i) ~= 0
        
        %change directory to path to look for folder
        cd(qfsmPackPaths{i,1})
        
        if actinConditions.maskLoc(i) == 1
            try
                % TN2022: if only External_ROI folder is provided
                masks = dir([qfsmPackPaths{i,1} fs char(extlroi(isfolder(extlroi)))]);
            catch ME
                switch ME.identifier
                    case 'MATLAB:catenate:dimensionMismatch'
                        % TN2023: Windows OS is not case-sensitive and will
                        % likely give multiple possible equivalent folder names.
                        tmpFolder = extlroi(isfolder(extlroi));
                        masks = dir([qfsmPackPaths{i,1} fs char(tmpFolder(1))]);
                    otherwise
                        % TN2019: these 2 folders contain essentially the same mask
                        % because to use the hand-drawn external mask, we have to pass
                        % that through SegmentationPackage.
                        masks = dir([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
                end
            end
        elseif actinConditions.maskLoc(i) == 2
            
            masks = dir([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
        elseif actinConditions.maskLoc(i) == 3
            
            masks = dir([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' fs 'refined_masks_for_channel_1']);
        else
            
            disp('error:: Folder of masks does not exist.')
        end
        
    else
        error('Please input cell mask.') % TN.AD. 20211119
    end
    %only retain movieInfo features that exist within mask
    
    %% Get speckle props
    %actinPropPerTimeInt{i,1}(1,1).speckleList = find(~isnan(MPM(:,2))); %%%% ISN'T THIS LINE IS USELESS? IF WE NEED TO INITIATE FIELD .SPECKLELIST,
    %%%% CAN'T WE JUST DO: actinPropPerTimeInt{i,1}(1,1).speckleList = [];
    
    for j = 1 : fsmFrames(i,1)
        
        %% get ILmax & IBkg - TN 20210524
        
        % load speckle cands structure for intensity of speckle and background
        if fsmFrames(i,1) < 10
            load([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'speckles' fs 'cands_' num2str(j,'%01.f') '.mat'], 'cands');
        elseif fsmFrames(i,1) < 100
            load([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'speckles' fs 'cands_' num2str(j,'%02.f') '.mat'], 'cands');
        elseif fsmFrames(i,1) < 1000
            % if actin movie is >100 frames then the first cands will be named ..001.mat and so on.
            load([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'speckles' fs 'cands_' num2str(j,'%03.f') '.mat'], 'cands');
        else
            error('Loading cands is not accounting for >999 fsm frames');
        end
        
        lmaxAll = vertcat(cands.Lmax);
                
        %% get .kinScorePos, .kinScore, .specklePosInit, .speckleVelocity, and
        %% .speckleSpeed
        
        %assigning speckle list, initial positions, velocities, and speeds
        actinPropPerTimeInt{i,1}(j,1).speckleList = find(~isnan(MPM(:,2*j)));
        
        %assigning positions of speckles with added noise (as (y,x) in
        %MPM's coordinate)
        actinPropPerTimeInt{i,1}(j,1).speckleInitPos = ...
            [MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j), ...
            MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j - 1)];
        
        if actinConditions.maskLoc(i) ~= 0
            
            if actinConditions.maskLoc(i) == 1
                
                              
                try
                    % these 2 folders (External_ROI and refined_masks from Segmentation Package contain essentially the same mask
                    
                    [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                        actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                        find(imread([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' ...
                        fs 'refined_masks_for_channel_1' fs masks(j + 2).name]));
                catch % if only External_ROI was provided
                    %extract (y,x) coordinates for all pixels within external mask
                    try
                        [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                            actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                            find(imread([qfsmPackPaths{i,1} fs char(extlroi(isfolder(extlroi))) fs ...
                            masks(j + 2).name]));
                    catch ME
                        switch ME.identifier
                            case 'MATLAB:catenate:dimensionMismatch'
                                tmpFolder = extlroi(isfolder(extlroi));
                                [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                                    actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                                    find(imread([qfsmPackPaths{i,1} fs char(tmpFolder(1)) fs ...
                                    masks(j + 2).name]));
                        end
                    end
                end
                
            elseif actinConditions.maskLoc(i) == 2
                
                %extract (y,x) coordinates for all pixels within refined mask
                [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                    actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                    find(imread([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'refined_masks' ...
                    fs 'refined_masks_for_channel_1' fs masks(j + 2).name]));
                
            elseif actinConditions.maskLoc(i) == 3
                
                [actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                    actinPropPerTimeInt{i,1}(j).indxMask(:,1)] = ...
                    find(imread([qfsmPackPaths{i,1} fs 'SegmentationPackage' fs 'refined_masks' ...
                    fs 'refined_masks_for_channel_1' fs masks(j + 2).name]));
            end
            
            %indices of positions
            pos = [actinPropPerTimeInt{i,1}(j).speckleInitPos(:,2) ...
                actinPropPerTimeInt{i,1}(j).speckleInitPos(:,1)];
            
            actinPropPerTimeInt{i,1}(j).speckleMaskDensity = ...
                sum(ismember(round(pos),[actinPropPerTimeInt{i,1}(j).indxMask(:,2), ...
                actinPropPerTimeInt{i,1}(j).indxMask(:,1)],'rows')) / ...
                length(actinPropPerTimeInt{i,1}(j).indxMask(:,1));
        end
        
        %vector of directories
        temp(j) = ...
            load(strcat([qfsmPackPaths{i,1} fs 'QFSMPackage' fs 'kineticAnalysis' fs maps{i}(j).name]));
        
        %assigning kinetic scores and their positions
        actinPropPerTimeInt{i,1}(j,1).kinScore = temp(j).kinScore(:,4);
        actinPropPerTimeInt{i,1}(j,1).kinScorePos = ...
            [temp(j).kinScore(:,3), temp(j).kinScore(:,2)];
        
        
        %assigning intensity of speckle and background
        ILmaxMat = nan(size(actinPropPerTimeInt{i,1}(j,1).speckleList,1),1);
        IBkgMat = nan(size(actinPropPerTimeInt{i,1}(j,1).speckleList,1),1);
            
        for iSpkl = 1:size(actinPropPerTimeInt{i,1}(j,1).speckleList,1) % loop through each speckle detected in each FSM frame
            
            xPos = actinPropPerTimeInt{i,1}(j,1).speckleInitPos(iSpkl,2);
            yPos = actinPropPerTimeInt{i,1}(j,1).speckleInitPos(iSpkl,1);
            indSpklCurr = find(lmaxAll(:,1) == xPos & lmaxAll(:,2) == yPos);
            
            % Check if speckle location from cands and MPM match (expect
            % everything to be matched, in the case that ghost speckles are removed)
            if isempty(indSpklCurr)
                error('Unknown error: Cands & MPM not matched. Speckle is not ghost!')
            end
            
            ILmax = cands(indSpklCurr).ILmax;
            IBkg  = cands(indSpklCurr).IBkg;
            
            % Save ILmax & IBkg in output structure
            ILmaxMat(iSpkl,1) = ILmax;
            IBkgMat(iSpkl,1) = IBkg;
            
        end
        actinPropPerTimeInt{i,1}(j,1).ILmax = ILmaxMat;
        actinPropPerTimeInt{i,1}(j,1).IBkg = IBkgMat;
        
        
        %calculate speckle velocities and speeds
        actinPropPerTimeInt{i,1}(j,1).speckleVelocity = ...
            [MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j + 2) - ...
            actinPropPerTimeInt{i,1}(j,1).speckleInitPos(:,1), ...
            MPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j + 1) - ...
            actinPropPerTimeInt{i,1}(j,1).speckleInitPos(:,2)];
        
        actinPropPerTimeInt{i,1}(j,1).speckleSpeed = ...
            sqrt(sum(actinPropPerTimeInt{i,1}(j,1).speckleVelocity .^ 2,2));
        
        %assigning speckle lifetime - TN20190828
        actinPropPerTimeInt{i,1}(j,1).speckleLifetime = ...
            lifetimeMPM(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j);
        
        %assigning speckle lifetime flag - KJ20210707
        actinPropPerTimeInt{i,1}(j,1).flagCompleteLft = ...
            flagCompleteLft(actinPropPerTimeInt{i,1}(j,1).speckleList,2*j);
        
        %assigning speckle midpoint positions
        actinPropPerTimeInt{i,1}(j,1).speckleMidPos = ...
            actinPropPerTimeInt{i,1}(j,1).speckleInitPos + ...
            (actinPropPerTimeInt{i,1}(j,1).speckleVelocity / 2);
        
        %distances of all speckles between each other for a given frame of speckles
        distMats = createDistanceMatrix(actinPropPerTimeInt{i}(j).speckleInitPos, ...
            actinPropPerTimeInt{i}(j).speckleInitPos);
        
        neighbhArea = (pi*(actinConditions.rad_density(i))^2);
        for iCol = 1 : size(distMats,2)
            matchindx = find(distMats(:,iCol) <= actinConditions.rad_density(i));
            
            actinPropPerTimeInt{i}(j).speckleDensity(iCol,1) = ...
                length(matchindx) / neighbhArea;
            
           matchedVect = actinPropPerTimeInt{i,1}(j,1).speckleVelocity(matchindx,:);
           matchedSpeed = actinPropPerTimeInt{i,1}(j,1).speckleSpeed(matchindx,:);
            
            %TN 20210808: add ILmaxNeighb & IBkgNeighb (to mirror
            %ilmaxCorrected and ibkgCorrected when 
            %    If no matched neighbor, ilmaxNeighb = ilmax;
            %    If there is a matched neighbor, ilmaxNeighb = (average ilmax)
            matchedILmax = actinPropPerTimeInt{i,1}(j,1).ILmax(matchindx,:);
            actinPropPerTimeInt{i}(j).ILmaxNeighb(iCol,1) = ...
                mean(matchedILmax);
            
            matchedIBkg = actinPropPerTimeInt{i,1}(j,1).IBkg(matchindx,:);
            actinPropPerTimeInt{i}(j).IBkgNeighb(iCol,1) = ...
                mean(matchedIBkg);
            
                        
            %KJ 20210712: if only one speckle neighbor exists, or only one
            %speckle neighbor has a velocity > 0, movement coherence is not
            %relevant
            if size(matchindx,1) == 1 || length(find(matchedSpeed > 0)) <= 1
                actinPropPerTimeInt{i}(j).speckleMvmtCohere(iCol,1) = NaN;
            else    
                
                % To calculate speckleMvmtCohere, take the resultant vector
                % magnitude, divide that by the average of speed of each
                % vector that was added - TN 20190913; sum-over-sum
                % strategy - TN 20191121
                
                resultantVect = sum(matchedVect);
                
                actinPropPerTimeInt{i}(j).speckleMvmtCohere(iCol,1) = ...
                    sqrt(sum(resultantVect.^2)) / sum(matchedSpeed);
                
            end
        end
        
    end %(for j = 1 : fsmFrames(i,1))
    
end %(for i = 1 : numMovs)

end

%% ~~~ the end ~~~
