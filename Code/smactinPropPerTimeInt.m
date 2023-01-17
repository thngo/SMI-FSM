function [combinePropPerTimeInt] = smactinPropPerTimeInt(smPropPerTimeInt,...
    actinPropPerTimeInt,combineConditions,threshMet,smactinFlag)
%SMACTINPROPPERTIMEINT : put together single molecule tracklets and associated local actin speckle properties
%
%SYNOPSIS: [combinePropPerTimeInt] = smactinPropPerTimeInt(smPropPerTimeInt,...
%    actinPropPerTimeInt,combineConditions,threshMet,smactinFlag)
%
%INPUT: smPropPerTimeInt  : (n x 1 cell array of structures) SMI arm
%                           analysis results of SMI-FSM analysis for n
%                           movies; contains various properties of SM. 
%                           (see output of smPropPerTimeInterval.m for details).
%
%     actinPropPerTimeInt : (n x 1 cell array of structures) FSM arm
%                           analysis results of SMI-FSM analysis for n
%                           movies; each structure contain speckle
%                           properties for m FSM intervals. 
%                           (see output of actinPropPerTimeInterval.m for details).
%
%       combineConditions : vector of structures with various fields imposing
%                           restrictions to associate SM and actin speckles.
%
%                       .matchRadius   : (number vector) radius of SM
%                                        neighborhood for matching to
%                                        neighboring speckles.
%                                       Example: .matchRadius = 3.5 (empirically 
%                                       determined from camera chromatic shift)
%
%                       .fsmInterval   : number of SM frames to be associated with 
%                                        FSM speckle frame (AKA per FSM frame).
%                                       Example: .fsmInterval = 50
%
%                       .firstFsmFrame : the first frame of speckles that should be
%                                        considered for matching to the start of SM movie.
%                                       Example: .firstFsmFrame = 2 (for
%                                       darkframe-scheme dataset)
%
%               threshMet : vector of cell vectors containing SM indices
%                           that survived a minimum lifetime (observations)
%                           threshold and exist within a mask obtained from
%                           qFSM.
%
%             smactinFlag : (structure) flags describing the method by
%                           which SM tracklets and actin speckles are
%                           associated spatially; contains the following
%                           field:
%
%                       .combFlag: integer with value either 7 (to match SM
%                                  to speckles based on average SM
%                                  position) or 10 (to match SM to speckles
%                                  based on frame-by-frame SM position).
%
%
%OUTPUT: combinePropPerTimeInt: Cell array with length = number of movies.
%                           Each cell contains two cells. First cell
%                           contains individual matched speckle
%                           information per single-molecule. Second cell
%                           contains average matched speckle information per
%                           single-molecule.
%
% NOTE: The following fields belong to structures contained within mentioned column
%       of cells; i.e. at end of each field description. Four categories of fields
%
% 1. Single-Molecules
%               smList: List of particles that exist in the window of time
%                       for time-averaging                          (c3,c6,c7,c8,c9)
%           .smMeanPos: Mean position (x/y-coordinates) per particle in
%                       the window for time-averaging               (c3,c6,c7,c8,c9)
%          .smMeanDisp: mean frame-to-frame displacement for all SM
%                       thresholding survivors in the window of
%                       time-averaging (units: pixels)              (c3,c6,c7,c8,c9)
%         .smNetVelocity: SM velocity vector  per frame averaged
%                       across a particular interval                (c3,c6,c7,c8,c9)
%          .smNetSpeed: net speed for all SM thresholding survivors
%                       within a particular interval.
%                       (units: pixels/SMframe)                     (c3,c6,c7,c8,c9)
%           .smMeanAmp: mean intensity of a SM track during a particular
%                       window of time                              (c3,c6,c7,c8,c9)
%          .trackClass: classification of all tracks based on moment scaling
%                       spectrum analysis
%                       = 0 : immobile
%                       = 1 : confined brownian
%                       = 2 : pure brownian (free diffusion)
%                       = 3 : directed motion                       (c3,c6,c7,c8,c9)
%            .diffCoef: diffusion coefficient with row number corresponding
%                       to particle index                           (c3,c6,c7,c8,c9)
%             .confRad: confinement radius of immobile, confined particle with row
%                       number corresponding to particle index      (c3,c6,c7,c8,c9)
%   .smComovementAngle: angle of comovement between the receptor and it nn
%                       speckle or sm                               (c3,c7,c8,c9)
%      .smIndiPosPerFr: list of SM position for each SM frame (group of ~51
%                       frames per interval)                        (c7)
%
% 2. Speckles
%         .speckleList: List of speckles that exist in the frame
%                                                                         (c1,c2,c3,c4,c7)
%          .specklePos: The speckle position for a particular FSM
%                       frame (units: pixels)                             (c1,c2,c3,c4,c7)
%     .speckleVelocity: The speckle velocity from a frame to the
%                       following frame (units: pixels / receptor frame)  (c1,c2,c3,c4,c7)
%        .speckleSpeed: speckle speed from a particular FSM frame to the
%                       following frame (units: pixels / receptor frame)  (c1,c2,c3,c4,c7)
%     .spklMtchPerSmFr: list of speckles matched per sm frame (group of ~51
%                       frames per interval)                              (c7)
%
% 3. Kinetic Scores
%         .kinScorePos: The positions at which we have kinetic scores
%                       (units: pixels)                             (c2,c4,c5,c6,c8)
%            .kinScore: The kinetic score from the qFSM package     (c2,c4,c5,c6,c8)
%
% 4. Miscellaneous
%         .maskDensity: The number of speckles within a mask divided
%                       by the number of pixels comprising the mask    (all columns)
%          .dist_match: Distance of a particular nn or collection of nns
%                       depending on which structure                   (all columns)
%                 .num: Number of neighbors contained within a specified
%                       disc                                           (all columns)
%           .rad_match: Matching radius to form disc and obtain neighbors
%                       for each speckle, kinetic score, or sm         (all columns)
%
%GLOSSARY:  SM  : single molecule
%           SMI : sinlge molecule imaging
%           QFSM: quantitative Fluorescent Speckle Microscopy
%           FL  : full-length (refer to compTrack pre-chopping)
%
% Deryl Tschoerner, February 2018
% Tra Ngo, modified Oct 2019
% Khuloud Jaqaman, July 2021, heavily modified to clean up and fix some bugs
%
%NOTE KJ 20210716: In the (near) future, change the output of
%smPropPerTimeInterval such that only the information of particles inside
%the mask is retained. When this change is made there, then this function
%must be changed so that it no longer needs/uses threshMet.
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



%% Const & variable initiation

numMovs = length(actinPropPerTimeInt);
fsmFrames = nan(numMovs,1);
smIntervals = nan(numMovs,1);

%% Output

combinePropPerTimeInt = cell(numMovs,1);

%% Calculation & matching:

for iMov = 1 : numMovs % loop through each movie
    
    % If no speckles or single-molecules in this movie: no analysis; this row (AKA cell) will be empty.
    if isempty(actinPropPerTimeInt{iMov,1}) || isempty(smPropPerTimeInt{iMov,1})
        continue
    end
    
    % intervals of frames per movie
    fsmFrames(iMov) = length(actinPropPerTimeInt{iMov});
    smIntervals(iMov) = length(smPropPerTimeInt{iMov});
    
    %% single molecule property assignment (+ speckle density within mask):
    
    %initialization
    combinePropPerTimeInt{iMov}{1} = repmat(struct('speckleMaskDensity',[],'rad_match', ...
        combineConditions.matchRadius(iMov),'num_match',[],'dist_match',[]),fsmFrames(iMov) - 1, 1);
    
    % Looping over the time windows to be considered.
    % Reminder: if combineConditions.firstFsmFrame = 2, then 1st SM
    % interval is skipped, since we pad the 1st interval to obtain
    % rate information for speckles (in previous analysis).
    for k = combineConditions.firstFsmFrame(iMov) : smIntervals(iMov)
        
        % density of speckles within mask
        combinePropPerTimeInt{iMov}{1}(k).speckleMaskDensity = ...
            actinPropPerTimeInt{iMov}(k).speckleMaskDensity;
        
        %single-molecule properties
        
        combinePropPerTimeInt{iMov}{1}(k).smList = ...
            threshMet{iMov}{k};
        combinePropPerTimeInt{iMov}{1}(k).smMeanPos = ...
            smPropPerTimeInt{iMov}(k).meanPos(threshMet{iMov}{k},:);
        combinePropPerTimeInt{iMov}{1}(k).smNetSpeed = ...
            smPropPerTimeInt{iMov}(k).netSpeed(threshMet{iMov}{k},:);
        combinePropPerTimeInt{iMov}{1}(k).smMeanAmp = ...
            smPropPerTimeInt{iMov}(k).meanAmp(threshMet{iMov}{k});
        combinePropPerTimeInt{iMov}{1}(k).smLifetime = ...
            smPropPerTimeInt{iMov}(k).lifetime(threshMet{iMov}{k});
        combinePropPerTimeInt{iMov}{1}(k).diffCoef = ...
            smPropPerTimeInt{iMov}(k).diffCoef(threshMet{iMov}{k});
        
        if isfield(smPropPerTimeInt{iMov}(k),'diffCoefF2F')
            combinePropPerTimeInt{iMov}{1}(k).smDiffCoefF2F = ...
                smPropPerTimeInt{iMov}(k).diffCoefF2F(threshMet{iMov}{k});
        else % .diffCoefF2F was created for Mode Analysis. If no mode analysis was implemented, make combinePropPerTimeInt.smDiffCoefF2F same as .diffCoef
            combinePropPerTimeInt{iMov}{1}(k).smDiffCoefF2F = ...
                combinePropPerTimeInt{iMov}{1}(k).diffCoef;
        end
        
        combinePropPerTimeInt{iMov}{1}(k).smDensity = ...
            smPropPerTimeInt{iMov}(k).smDensity(threshMet{iMov}{k});
        combinePropPerTimeInt{iMov}{1}(k).smIndiPosPerFr = ...
            smPropPerTimeInt{iMov}(k).indiPosPerFr(threshMet{iMov}{k});
        
        if isfield(smPropPerTimeInt{iMov}(k),'diffRad')
            combinePropPerTimeInt{iMov}{1}(k).smDiffRad = ...
                smPropPerTimeInt{iMov}(k).diffRad(threshMet{iMov}{k});
        else
            combinePropPerTimeInt{iMov}{1}(k).smDiffRad = ...
                smPropPerTimeInt{iMov}(k).confRad(threshMet{iMov}{k}); % make diffRad similar to confRad
        end
        
        if isfield(smPropPerTimeInt{iMov}(k),'diffMode')
            combinePropPerTimeInt{iMov}{1}(k).smDiffMode = ...
                smPropPerTimeInt{iMov}(k).diffMode(threshMet{iMov}{k});
        end
        
        if isfield(smPropPerTimeInt{iMov}(k),'aggregState') % check if .aggregState is a field for backward compartibility b/c it's only added later on
            combinePropPerTimeInt{iMov}{1}(k).smAggregState = ...
                smPropPerTimeInt{iMov}(k).aggregState(threshMet{iMov}{k});
        end
        
        combinePropPerTimeInt{iMov}{1}(k).smMeanDisp = ...
            smPropPerTimeInt{iMov}(k).meanDisp(threshMet{iMov}{k});
        combinePropPerTimeInt{iMov}{1}(k).smNetVelocity = ...
            smPropPerTimeInt{iMov}(k).netVelocity(threshMet{iMov}{k},:);
        combinePropPerTimeInt{iMov}{1}(k).confRad = ...
            smPropPerTimeInt{iMov}(k).confRad(threshMet{iMov}{k});
        combinePropPerTimeInt{iMov}{1}(k).trackClass = ...
            smPropPerTimeInt{iMov}(k).trackClass(threshMet{iMov}{k});
        
    end %(for k = combineConditions.firstFsmFrame(iMov) : smIntervals(iMov))
    
    combinePropPerTimeInt{iMov}{2} = combinePropPerTimeInt{iMov}{1};
    
    %% nn speckle property assignment:
    
    for k = combineConditions.firstFsmFrame(iMov): smIntervals(iMov) % looping over time windows
        
        smMatchPt = smPropPerTimeInt{iMov}(k).meanPos(threshMet{iMov}{k},:);
        actinMatchPt = actinPropPerTimeInt{iMov}(k).speckleInitPos;
        
        % if no speckle or no sm in this time window, then no analysis and go to the next time window.
        if isempty(actinMatchPt) || isempty(smMatchPt)
            continue
        end
        
        distMat = createDistanceMatrix(actinMatchPt, smMatchPt);
        
        for iCol = 1 : size(distMat,2) % given SM objects, find speckle nns and get their properties
            
            switch smactinFlag.combFlag(iMov)
                
                case 7 % match neighboring speckles based on average SM position
                    
                    nn = [];
                    mtchindx = find(distMat(:,iCol) <= combinePropPerTimeInt{iMov}{1}(k).rad_match);
                    [mtchIndxBySmFr,matSpkAtSmFrIndx] = deal(mtchindx); %KJ 20210709: Add this here, because otherwise below code will crash for this case
                    
                case 10 % match neighboring speckles based on frame-by-frame SM position
                    
                    nn = [];
                                        
                    trackDistMat = createDistanceMatrix(actinPropPerTimeInt{iMov}(k).speckleInitPos, ...
                        combinePropPerTimeInt{iMov}{1}(k).smIndiPosPerFr{iCol,1});
                    mtchSmFrIndxMat = trackDistMat <= combinePropPerTimeInt{iMov}{1}(k).rad_match;
                    
                    maxNumSpek = max(sum(mtchSmFrIndxMat)); % get max number of speckles matched at a single frame
                    numSmFr = size(trackDistMat, 2); % get SM track lifetime
                    
                    matSpkAtSmFr = nan(maxNumSpek, numSmFr); % initialize matrix containing matched speckles
                    
                    for iTime = 1:numSmFr
                        spekInd = find(mtchSmFrIndxMat(:,iTime) == 1);
                        matSpkAtSmFr(1:length(spekInd),iTime) = ...
                            actinPropPerTimeInt{iMov}(k).speckleList(spekInd);
                    end
                    
                    combinePropPerTimeInt{iMov}{1}(k).([nn 'spklMtchPerSmFr']){iCol,1} = matSpkAtSmFr;
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'spklMtchPerSmFr']){iCol,1} = matSpkAtSmFr;
                    
                    mtchList = unique(matSpkAtSmFr);
                    mtchListLogic = ~isnan(mtchList);
                    mtchID = mtchList(mtchListLogic);
                    mtchID = mtchID(:); %KJ 200204: make sure that mtchID is always a column vector
                    
                    mtchindx = nan(size(mtchID));
                    for iMtchIndx = 1:length(mtchID)
                        mtchindx(iMtchIndx) = ...
                            find(actinPropPerTimeInt{iMov}(k).speckleList == mtchID(iMtchIndx));
                    end
                    
                    %KJ 20210712: get speckle indices first in matrix
                    %format of speckle matching (per frame)
                    matSpkAtSmFrIndx = NaN(size(matSpkAtSmFr));
                    indxGood = find(~isnan(matSpkAtSmFr));
                    indxGood = indxGood(:)';
                    for iSpk = indxGood
                        matSpkAtSmFrIndx(iSpk) = find(actinPropPerTimeInt{iMov}(k).speckleList == matSpkAtSmFr(iSpk));
                    end
                    
                    % TN 20210601: get full list of speckles matched to sm by frame (with repeat, if
                    % sm is matched to the same speckle multiple times)
                    %KJ 20210712: Because mapping has already been done at
                    %matrix level, no need for mapping here now

                    mtchIndxBySmFr = matSpkAtSmFrIndx(~isnan(matSpkAtSmFrIndx));
                    mtchIndxBySmFr = mtchIndxBySmFr(:);
                   
            end
            
            % Retrieve indices of nearest neighbors
            switch smactinFlag.combFlag(iMov)
                
                case 7 % assign nn speckle indices based on average SM position
                    combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleList']){iCol,1} = ...
                        actinPropPerTimeInt{iMov}(k).speckleList(mtchindx);
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleList']){iCol,1} = ...
                        combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleList']){iCol,1};
                    
                case 10 % assign nn speckle indices based on frame-by-frame matching
                    combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleList']){iCol,1} = ...
                        mtchID;
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleList']){iCol,1} = ...
                        mtchID;
            end
            
            % Distances of neighbors:
            combinePropPerTimeInt{iMov}{1}(k).dist_match{iCol,1} = distMat(mtchindx,iCol); %KJ 200204 comment: This is distance from mean position, even when matching is done based on frame-by-frame position
            
            % Number of neighbors
            combinePropPerTimeInt{iMov}{1}(k).num_match(iCol,1) = ...
                length(combinePropPerTimeInt{iMov}{1}(k).dist_match{iCol,1});
            
            % Assignment of neighboring speckle properties to the SM object
            % Changed average to weighted average (Tra 20210602 and KJ 20210708-20210712)
            
            %Position
            
            combinePropPerTimeInt{iMov}{1}(k).([nn 'specklePos']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).speckleInitPos(mtchindx,:);
            combinePropPerTimeInt{iMov}{2}(k).([nn 'specklePos'])(iCol,:) = ...
                [mean(actinPropPerTimeInt{iMov}(k).speckleInitPos(mtchIndxBySmFr,1)), ...
                mean(actinPropPerTimeInt{iMov}(k).speckleInitPos(mtchIndxBySmFr,2))];
            
            %Movement
            
            combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleVelocity']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).speckleVelocity(mtchindx,:);
            combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleVelocity'])(iCol,:) = ...
                [nanmean(actinPropPerTimeInt{iMov}(k).speckleVelocity(mtchIndxBySmFr,1)), ...
                nanmean(actinPropPerTimeInt{iMov}(k).speckleVelocity(mtchIndxBySmFr,2))];
            
            combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleSpeed']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).speckleSpeed(mtchindx) ...
                / combineConditions.fsmInterval(iMov);
            combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleSpeed'])(iCol,1) = ...
                sqrt(sum(combinePropPerTimeInt{iMov}{2}(k).speckleVelocity(iCol,:) .^ 2)) ...
                / combineConditions.fsmInterval(iMov); % net speed
            
            %KJ 20210712: the above speed is net speed, i.e. magnitude of
            %net velocity. This here is the average speed regardless of
            %direction. This field is needed to calculate movement
            %coherenece parameter below.
            speckleSpeedVec = actinPropPerTimeInt{iMov}(k).speckleSpeed(mtchIndxBySmFr,1);
            combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleSpeedAveMag'])(iCol,1) = ...
                nanmean(speckleSpeedVec) / combineConditions.fsmInterval(iMov);

            combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleMvmtCohere']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).speckleMvmtCohere(mtchindx);
            
            %KJ 20210709: if there is only one speckle neighbor, then
            %speckle movement coherence is not relevant
            %20210712: same if there is at most one speckle with nonzero velocity
            if length(mtchindx) == 1 || length(find(speckleSpeedVec>0)) <= 1 
                combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleMvmtCohere'])(iCol,1) = NaN;
            else
                combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleMvmtCohere'])(iCol,1) = ...
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleSpeed'])(iCol,1) ...
                    / combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleSpeedAveMag'])(iCol,1);
            end
            
            % KJ 20210709: Angle between speckle velocity (individual or
            % net) and SM velocity - Major bug fixes
            if isempty(mtchindx)
                
                combinePropPerTimeInt{iMov}{1}(k).([nn 'smComvmtAngle']){iCol,1} = [];
                combinePropPerTimeInt{iMov}{2}(k).([nn 'smComvmtAngle'])(iCol,:) = NaN;
                
            else
                
                vel1 = combinePropPerTimeInt{iMov}{1}(k).smNetVelocity(iCol,:);
                normVel1 = norm(vel1);
                for iAng = 1 : length(mtchindx)
                    vel2 = combinePropPerTimeInt{iMov}{1}(k).speckleVelocity{iCol,1}(iAng,:);
                    normVel2 = norm(vel2);
                    if normVel2 > 0 && normVel1 > 0
                        combinePropPerTimeInt{iMov}{1}(k).([nn 'smComvmtAngle']){iCol,1}(iAng,1) = ...
                            acosd(dot(vel1,vel2) / (normVel1*normVel2));
                    else
                        combinePropPerTimeInt{iMov}{1}(k).([nn 'smComvmtAngle']){iCol,1}(iAng,1) = NaN;
                    end
                end
                
                vel2 = combinePropPerTimeInt{iMov}{2}(k).speckleVelocity(iCol,:);
                normVel2 = norm(vel2);
                if normVel2 > 0 && normVel1 > 0
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'smComvmtAngle'])(iCol,1) = ...
                        acosd(dot(vel1,vel2) / (normVel1*normVel2));
                else
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'smComvmtAngle'])(iCol,1) = NaN;
                end
                
            end
            
            %Intensity
            
            %KJ 20210712: Fixed bug in calculating ILMaxCorrected and
            %IBkgCorrected. Previous implementation (in 20210602) did not
            %handle properly the case when multiple speckles are matched to
            %an SM in a single frame

            intPerFrame = NaN(size(matSpkAtSmFrIndx));
            intPerFrame(indxGood) = actinPropPerTimeInt{iMov}(k).ILmax(matSpkAtSmFrIndx(indxGood)); %indxGood is from Line 276 above
            intPerFrameMean = nanmean(intPerFrame,1);
            
            combinePropPerTimeInt{iMov}{1}(k).([nn 'ILmax']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).ILmax(mtchindx);            
            combinePropPerTimeInt{iMov}{2}(k).([nn 'ILmax'])(iCol,1) = ...
                nanmean(intPerFrameMean);
            switch smactinFlag.combFlag(iMov)                
                case 7
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'ILmaxCorrected'])(iCol,1) = ...
                        combinePropPerTimeInt{iMov}{2}(k).([nn 'ILmax'])(iCol,1);
                case 10
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'ILmaxCorrected'])(iCol,1) = ...
                        nansum(intPerFrameMean)/numSmFr;
            end
            
            intPerFrame = NaN(size(matSpkAtSmFrIndx));
            intPerFrame(indxGood) = actinPropPerTimeInt{iMov}(k).IBkg(matSpkAtSmFrIndx(indxGood));
            intPerFrameMean = nanmean(intPerFrame,1);

            combinePropPerTimeInt{iMov}{1}(k).([nn 'IBkg']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).IBkg(mtchindx);
            combinePropPerTimeInt{iMov}{2}(k).([nn 'IBkg'])(iCol,1) = ...
                nanmean(intPerFrameMean);
            switch smactinFlag.combFlag(iMov)
                case 7
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'IBkgCorrected'])(iCol,1) = ...
                        combinePropPerTimeInt{iMov}{2}(k).([nn 'IBkg'])(iCol,1);
                case 10
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'IBkgCorrected'])(iCol,1) = ...
                        nansum(intPerFrameMean)/numSmFr;
            end
            
            %Density
            
            combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleDensity']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).speckleDensity(mtchindx);
            
            switch smactinFlag.combFlag(iMov)
                
                case 7  
                    
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleDensity'])(iCol,1) = ...
                        size(mtchindx,1)/(pi*(combinePropPerTimeInt{iMov}{2}(k).rad_match)^2);
                    
                case 10
                    
                    spklCountEachFr = nansum(~isnan(combinePropPerTimeInt{iMov}{1}(k).spklMtchPerSmFr{iCol,1}),1);
                    meanSpklCount = mean(spklCountEachFr);
                    searchArea = pi*(combinePropPerTimeInt{iMov}{2}(k).rad_match)^2;
                    
                    combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleDensity'])(iCol,1) = ...
                        meanSpklCount/searchArea;
                    
            end
            
            %Lifetime
            
            combinePropPerTimeInt{iMov}{1}(k).([nn 'speckleLifetime']){iCol,1} = ...
                actinPropPerTimeInt{iMov}(k).speckleLifetime(mtchindx);
            combinePropPerTimeInt{iMov}{2}(k).([nn 'speckleLifetime'])(iCol,1) = ...
                nanmean(actinPropPerTimeInt{iMov}(k).speckleLifetime(mtchIndxBySmFr));
                        
        end %(for iCol = 1 : size(distMat,2))
        
    end %(for k = combineConditions.firstFsmFrame(iMov): smIntervals(iMov))
    
end %(for iMov = 1 : numMovs)

end
