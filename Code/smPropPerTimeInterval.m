function [smPropPerTimeInt,threshMet, summaryNumTrackAll, numPassDecomptrackAll] = ...
    smPropPerTimeInterval(utrackPackPaths, smConditions,qfsmPackPaths, ...
    diffModeDividerStruct, smSpanRadius, intensityInfoBackUp)
%SMPROPPERTIMEINTERVAL reads data from u-track analysis folders and calculates spatiotemporal properties of single-molecule data
%
%SYNOPSIS: [smPropPerTimeInt,threshMet, summaryNumTrackAll, numPassDecomptrackAll] = smPropPerTimeInterval(utrackPackPaths, smConditions, qfsmPackPaths, diffModeDividerStruct, smSpanRadius, intensityInfoBackUp)
%
%INPUT:  utrackPackPaths: (n x 1 cell) paths leading to n U-Track analysis
%                         results in "TrackingPackage" folders.
%
%           smConditions: (structure) parameters for processing tracks;
%                         contains the following fields.
%
%                       .minLft: (n x 1 vector of integer) minimum acceptable
%                                lifetime for single-molecules for motion 
%                                analysis. Any SM with lifetime lower than 
%                                the minLft will be eliminated from processing.
%                                Recommended value: 10 
%
%                       .rad_density: (n x 1 vector of double) radius used
%                                for calculating local SM density.
%                                Recommended value: 7 (empirically determined by Deryl T. 2018)
%
%                       .windAvg: (n x 1 vector of integer) number of SMI frames
%                                 used to average (split) SMI data into FSM intervals.
%
%                       .divFlag: (n x 1 vector of integer) flags
%                                 indicating how SMI frames are matched
%                                 temporally to FSM frames.
%                          = -1 : link SM frames to earlier FSM frames
%                          =  0 : link SM frames to center FSM frames
%                          =  1 : link SM frames to later FSM frames
%
%                       .maskLoc: (n x 1 vector of integer) flag indicating
%                                location of cell ROI mask.
%                          = 0 : there do not exist any masks
%                          = 1 : there exists a folder with external
%                                masks within path provided in
%                                qfsmPackPaths (make sure folder name
%                                follows a name in extlroi variable.)
%                          = 2 : there exists a folder with refined
%                                masks within the QFSMPackage folder
%                                provided in qfsmPackPaths
%                          = 3 : there exists a folder with refined
%                                masks within the SegmentationPackage
%                                folder provided in qfsmPackPaths
%                          if .maskLoc = [], the code will crash.
%
%          qfsmPackPaths: (n x 1 cell) paths leading to QFSM Package results
%
%           smSpanRadius: (integer) radius of domain that a single sm can
%                         span before being chopped into (multiple)
%                         spatially localized tracklets, with equal lifetimes.
%                         Optional. Default = 300/9
%
%  diffModeDividerStruct: (structure -- input not implemented) parameters
%                         for performing diffusion mode analysis; contains
%                         the following fields:
%                      .trajLength: Trajectory lengths for which mode divider values
%                        are supplied.
%                      .coordStd  : Coordinate standard deviations for which mode
%                        divider values are supplied.
%                      .divider   : (number of trajectory lengths) x (number of
%                        coordinates stds) x (number of modes - 1) array
%                        of mode divider values.
%                       Optional. Input not implemeted. Default = [].
%
%           intensityInfoBackUp: (optional, row vector of 2 values) data of
%                         monomer intensity as [dataMean dataStdev], used to
%                         estimate oligomeric state if automatic estimation
%                         at each movie error out. This option should only
%                         only be used to test movie with only 1 track,
%                         which did not contain enough info for estimating
%                         monomer intensities. -- TN 20210408
%
%OUTPUT:   smPropPerTimeInt: (n x 1 cell array of structures) SMI arm
%                            analysis results of SMI-FSM analysis for n
%                            movies; contains various properties of SM;
%                            each structure contains the following fields: 
%
%                   .smFrames        : List of SM frames associated with each FSM frame.
%                                  
%                   .smList          : List of SM tracks that exist in each FSM interval.
%
%                   .meanPos         : Mean position (x/y-coordinates) per SM track each FSM interval.
%
%                   .meanDisp        : Mean frame-to-frame displacement per SM track each FSM interval.
%
%                   .lifetime        : final frame of the last observation for a SM track
%                                  minus the first frame of the first observation for
%                                  a SM track within a particular FSM interval.
%
%                   .numObservations : number of observations for a SM track within a FSM interval.
%
%                   .netVelocity     : SM velocity vector  per frame averaged
%                                  across each FSM interval.
%
%                   .meanAmp         : mean intensity of a SM track during an FSM interval.
%
%                   .netSpeed        : the net speed for a particular SM track during an FSM interval.
%
%                   .netVelAngle     : the net velocity angle within an FSM interval.
%                                  For a SM track, this is the angle between .netSpeed 
%                                  direction and a horizontal axis crossing the initial
%                                  position of the SM track.
%
%                   .trackClass      : classification of all tracks based on moment 
%                                  scaling spectrum analysis.
%                                   = 0 : immobile
%                                   = 1 : confined brownian
%                                   = 2 : pure brownian (free diffusion)
%                                   = 3 : directed motion
%
%                   .diffCoef        : diffusion coefficient of each SM track in an FSM interval; 
%                                  row corresponding to SM index.
%
%                   .confRad         : confinement radius of immobile & confined SM track; 
%                                  row corresponds to SM index.
%
%                   .diffMode        : SM diffusion mode type as analysed by diffusionModeAnalysis
%
%                   .diffCoefF2F     : SM diffusion coefficient calculated from frame-to-frame displacement
%
%                   .diffRad         : SM diffusion radius from diffusionModeAnalysis (not implemented TN20210310)
%
%
%                   threshMet  : (n x 1 cell vectors) contains SM indices that survive 
%                             a minimum lifetime (observations) and exist within 
%                             cell ROI masks from qFSM.
%
%               summaryNumTrackAll: (table of n*m rows) summary of the number of 
%                                   tracklets added by each SM track processing 
%                                   steps for m FSM intervals of n movies;
%                                   contains the following columns: 
%                                           iMov
%                                           iFrame (indicates the FSM interval)
%                                           numTrack
%                                           trackWithTransMot
%                                           totTransMotTracklet
%                                           addedByDcmss
%                                           numTrackDcmss
%                                           addedByTime
%                                           numAfterTime
%                                           addedBySpan
%                                           numAfterSpan
%                                           numBigSpan
%
%            numPassDecomptrackAll: (table of n*m rows) summary of number
%                                   of tracks from u-track that is split by
%                                   each FSM interval; contains the 
%                                   following columns (added 20220714)
%                                           iMov
%                                           sizeMinDecompTrack
%                                           numPassTrack
%
%GLOSSARY:  SM  : single molecule
%           SMI : sinlge molecule imaging
%           QFSM: quantitative Fluorescent Speckle Microscopy
%           FL  : full-length (refer to compTrack pre-chopping)
%
% Deryl Tschoerner, February 2018
% Modified by Khuloud Jaqaman, February 2019
% Extensively modified by Tra H. Ngo, September 2019
%
%NOTE KJ 20210716: In the (near) future, change the output of this code
%such that only the information of particles inside the mask is retained.
%Then, there is no need for the output variable threshMet. When this change
%is made here, then smactinPropPerTimeInt must be changed so that it no
%longer needs/uses threshMet.
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
if isempty(utrackPackPaths)
    error('error:: single-particle tracking path was not found')
elseif isempty(qfsmPackPaths)
    error('error:: speckle information path was not found')
end

if nargin < 3
    error('error:: invalid number of input')
elseif nargin == 3
    diffModeDividerStruct = []; % if no input for diffModeDividerStruct, make this empty.
end

if ~exist('smSpanRadius','var')|| isempty(smSpanRadius)
    smSpanRadius = 300/9; % (pixels ~ 3um for pixel size = 90nm)
end

if ~exist('intensityInfoBackUp','var') || isempty(intensityInfoBackUp)
    intensityInfoBackUp = [0.0014,0.0004]; % some value from an arbitrary CD36 movie
end

%% Initialization & definitions
numMovs = length(utrackPackPaths); % number of SM movies being processed
fs = filesep; % define file separator (for different operating systems)
numFramesFSM = nan(numMovs,1); % number of FSM speckle frame
numFramesSM = nan(numMovs,1); % number of SM frame
extlroi = ["external_roi","external_mask1","External_mask","External_ROI","external_ROI", ...
    "External_Segmentation"]; % vector containing variations in naming of external mask folders

%% Output:
smPropPerTimeInt = cell(numMovs,1); % cell array containing SM physical metrics to be calculated within this function
threshMet = cell(numMovs,1); % cell array containing SM indices meeting criteria set by smCondition input
summaryNumTrackAll = []; % initialize table
numPassDecomptrackAll = [];  % initialize table
%% Looping through each sm movie
for iMov = 1 : numMovs
    
    % if one of the paths are empty, move on to the next listed path
    if isempty(utrackPackPaths{iMov,1})
        continue
    end
    
    % define smPropPerTimeInt output: each cell contains SM properties calculated for each movie
    smPropPerTimeInt{iMov,1} = struct;
    
    %% load SM tracks and relevant mask information
    
    % load qFSM movie mask (name, folder, etc) to exclude irrelevant SM
    if ~isempty(qfsmPackPaths{iMov,1}) && smConditions.maskLoc(iMov) ~= 0
        cd(qfsmPackPaths{iMov,1}) % change directories to search for the mask
        
        if smConditions.maskLoc(iMov) == 1
            try % works with Linux (case-sensitive)
                masks = dir([qfsmPackPaths{iMov,1} fs char(extlroi(isfolder(extlroi)))]);
            catch % works with Windows (in non-case-sensitive settings)
                tmpFold = dir(fullfile(qfsmPackPaths{iMov,1}, 'ext*'));
                masks = dir(fullfile(qfsmPackPaths{iMov,1}, tmpFold.name));
            end
        elseif smConditions.maskLoc(iMov) == 2
            masks = dir([qfsmPackPaths{iMov,1} fs 'QFSMPackage' fs ...
                'refined_masks' fs 'refined_masks_for_channel_1']);
        elseif smConditions.maskLoc(iMov) == 3
            masks = dir([qfsmPackPaths{iMov,1} fs 'SegmentationPackage' fs ...
                'refined_masks' fs 'refined_masks_for_channel_1']);
        else
            error('Folder of masks does not exist.')
        end
        
    end
    
    % load uTrack output in "TrackingPackage/tracks" folder
    try
        load([utrackPackPaths{iMov,1} fs 'tracks' fs ...
            'Channel_1_tracking_result.mat'], 'tracksFinal');
    catch
        error('tracksFinal DOES NOT LOAD!!!!!!!!')
    end
    
    %% Eliminating obviously useless data
    
    criteria.lifeTime.min = smConditions.minLft(iMov);
    indxKeep = chooseTracks(tracksFinal,criteria);
    tracksFinal = tracksFinal(indxKeep);
    
    %% Calculate full-length (FL) SM data size
    
    %KJ 20210714: get number of SM frames from sequence of events
    seqOfEvents = vertcat(tracksFinal.seqOfEvents);
    numFramesSM(iMov) = max(seqOfEvents(:,1));
    
    % calculate relevant interval between FSM and SM movies based on smConditions input
    numFramesFSM(iMov) = floor(numFramesSM(iMov) / smConditions.windAvg(iMov) + 1);
    
    %% Time frame assignment:: getting .smFrames
    
    % assigning which SM frame(s) belong to each FSM frame and which particle(s)
    % belong to which interval of frames - depending on divFlag input
    
    switch smConditions.divFlag(iMov)
        
        case 0
            
            %to divide an interval for centering FSM frame
            cHalf = ceil(smConditions.windAvg(iMov) / 2);
            
            %the first interval will begin with the first FSM frame
            windFirst = 1;
            
            %one more interval for centered
            windAdd = 1;
            
            %all intervals but first
            for iFrame = numFramesFSM(iMov) : -1 : 2
                smPropPerTimeInt{iMov,1}(iFrame,1).smFrames = ...
                    ((iFrame - 2)*smConditions.windAvg(iMov) + cHalf + 1: (iFrame - 1)*smConditions.windAvg(iMov) + cHalf + 1)';
            end
            
            %first interval
            smPropPerTimeInt{iMov,1}(1,1).smFrames = (1 : cHalf + 1)';
            
            %last interval
            indxKeep = find(smPropPerTimeInt{iMov,1}(numFramesFSM(iMov),1).smFrames <= numFramesSM(iMov));
            smPropPerTimeInt{iMov,1}(numFramesFSM(iMov),1).smFrames = smPropPerTimeInt{iMov,1}(numFramesFSM(iMov),1).smFrames(indxKeep); %#ok<FNDSB>
            
            
        case -1
            
            %no interval increment if not centered
            windAdd = 0;
            
            %the first interval will begin with the first FSM frame
            windFirst = 1;
            
            %memory allocation for frame intervals via assignment of last
            %interval
            lastFrameSM = min( (numFramesFSM(iMov)-1)*smConditions.windAvg(iMov)+1 , numFramesSM(iMov) );
            smPropPerTimeInt{iMov,1}(windFirst + (numFramesFSM(iMov) - 2),1).smFrames = ...
                ((numFramesFSM(iMov) - 2)*smConditions.windAvg(iMov) + 1 : lastFrameSM)';
            
            %assigning all other intervals
            for iFrame = windFirst : windFirst + (numFramesFSM(iMov) - 3)
                smPropPerTimeInt{iMov,1}(iFrame,1).smFrames = ...
                    ((iFrame - windFirst)*smConditions.windAvg(iMov) + 1 : ...
                    (iFrame + 1 - windFirst)*smConditions.windAvg(iMov) + 1)';
            end
            
        case 1
            
            %no interval increment if not centered
            windAdd = 0;
            
            %the first interval will begin with the second FSM frame
            windFirst = 2;
            
            %memory allocation for frame intervals via assignment of last
            %interval
            lastFrameSM = min( (numFramesFSM(iMov)-1)*smConditions.windAvg(iMov)+1 , numFramesSM(iMov) );
            smPropPerTimeInt{iMov,1}(windFirst + (numFramesFSM(iMov) - 2),1).smFrames = ...
                ((numFramesFSM(iMov) - 2)*smConditions.windAvg(iMov) + 1 : lastFrameSM)';
            
            %assigning all other intervals
            for iFrame = windFirst : windFirst + (numFramesFSM(iMov) - 3)
                smPropPerTimeInt{iMov,1}(iFrame,1).smFrames = ...
                    ((iFrame - windFirst)*smConditions.windAvg(iMov) + 1 : ...
                    (iFrame + 1 - windFirst)*smConditions.windAvg(iMov) + 1)';
            end
            
    end
    
    %% Clean up tracks
    
    [tracksReform,~] = removeSimultaneousMergeSplit(tracksFinal);
    [tracksReformat] = reformatSplitsMergesKeepLongerLivedSegment(tracksReform);
    thresholdTime = [];
    tracksClean = removeSplitMergeArtifactsChronological(tracksReformat,thresholdTime); % use default for 2nd input
    
    %% Assign clustering/oligomerization/assembly state
    
    try
        
        % Get intensity information (Monomer intensity mean & stdev) for each movie
        tmpTrackSEL = getTrackSEL(tracksClean);
        tmpRealStart = min(tmpTrackSEL (:,1)); % get earliest track start
        intensityFitEndFr = tmpRealStart + 9;
        intensityFitBegFr = tmpRealStart + 5;
        observations = getIntensityByFrameFromTrack(tracksClean, intensityFitEndFr, intensityFitBegFr); % collect observations (vector of intensity from frame range of tracksFinal)
        
        alpha = 1; % (use BIC instead of F-test)
        variableMean = -2;
        variableStd = 3; % by Mutch et al. Biophys 2007 (mean and standard deviation of mode n were ~ n times those of mode 1)
        showPlot = 0;
        numModeMinMax = 4; % empirically decided; (most CD36 fits don't even get to tetramer, increasing the maximum number of possible fitted mode doesn't change the fit result)
        binStrategy = 2; % default
        plotName = 'Figure';
        logData = 1; % log-normal fit (instead of Gaussian fit)
        modeParamIn = []; % default
        ratioTol = 0.1; % (empirically determined after looking at CD36 & VEGFR2 results)
        
        [~,~,modeParam,~,~] = ...
            fitHistWithGaussians(observations,alpha,variableMean,variableStd,...
            showPlot,numModeMinMax,binStrategy,plotName,logData,modeParamIn,ratioTol);
        
        if size(modeParam,1) > 1 % take mode 2 as the reference (mode 1 = 2nd mode's intensity divided by 2)
            
            param = modeParam(2,1:2); % log normal distribution parameters of mode 2
            dataMean = exp(param(1)+param(2)^2/2); % convert to data parameters
            %dataVar = exp(param(2)^2+2*param(1))*(exp(param(2)^2)-1);
            
            paramMode1 = modeParam(1,1:2);
            dataVarMode1 = exp(paramMode1(2)^2+2*paramMode1(1))*(exp(paramMode1(2)^2)-1);
            dataStdev1 = sqrt(dataVarMode1);
            
            intensityInfo = [dataMean/2 dataStdev1];
            % mean intensity of mode 1 = mean intensity of mode 2 / 2
            % we are taking the variance of mode 1 as fitted, because right now
            % this input is inconsequent for aggregStateFromCompTracksMIQP.
            
        else % Exception handling if the fit only has Mode 1
            
            param = modeParam(1,1:2); % log normal distribution parameters of mode 1
            dataMean = exp(param(1)+param(2)^2/2); % convert to data parameters
            dataVar = exp(param(2)^2+2*param(1))*(exp(param(2)^2)-1);
            dataStdev = sqrt(dataVar);
            
            intensityInfo = [dataMean dataStdev]; % monomer's mean
        end
        
    catch
        intensityInfo = intensityInfoBackUp;
        warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        warning('Intensity analysis error. Using back-up intensity information.')
        warning('This should only be allowed if user do not care about intensity information')
        warning('This should only be allowed if there is only 1 compound track from u-track.')
        warning('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
    
    % Get compTrack with oligomeric states
    [compTracksOut] = aggregStateFromCompTracksMIQP(tracksClean,intensityInfo);
    compTrackAltAgg = compTracksOut.alternativeFormatTracks;
    clear tracksClean
    
    % Decompounded-tracks
    [decompTrack] = decompoundCompTracks(compTrackAltAgg);
    
    % Count number of tracks that will get split by FSM intervals - TN20220714
    try
        tmpSmFr = [smPropPerTimeInt{iMov,1}.smFrames]; % matrix of all SMI frames, col = FSM intervals
        fsmSplitFrs = tmpSmFr(end,1:end-1);
    catch
        fsmSplitFrs = [51   101   151   201   251   301   351   401   451   501   551   601   651   701   751   801   851];
    end
    [~, numPassTrack] = selTrackChecker(decompTrack, fsmSplitFrs);
    
    [~, sizeMinDecompTrack] = selTrackChecker(decompTrack, [], smConditions.minLft(iMov));
    cntDecompTrack = table(iMov, sizeMinDecompTrack, numPassTrack);
    numPassDecomptrackAll = vertcat(numPassDecomptrackAll, cntDecompTrack);
    
    %% Loop over FSM intervals
    
    windLast = windFirst + windAdd + numFramesFSM(iMov) - 2; % last FSM interval
    
    for iFrame = windFirst : windLast % looping through each grouped FSM interval (e.g. 50 SM frame are grouped into 1 qFSM frame)
        
        %% Mask coordinates
        
        % get mask for each interval from mask of each corresponding FSM frame
        maskCoord = cell(windLast - windFirst + 1,1);
        
        % extract (y,x) coordinates for all pixels within refined, external, or segmented mask
        if smConditions.maskLoc(iMov) ~= 0
            if smConditions.maskLoc(iMov) == 1
                try
                    imgMask = imread([qfsmPackPaths{iMov,1} fs char(extlroi(isfolder(extlroi))) fs ...
                        masks(iFrame + 2).name]);
                    [maskCoord{iFrame,1}(:,2),maskCoord{iFrame,1}(:,1)] = ...
                        find(imgMask);
                catch
                    imgMask = imread(fullfile(masks(iFrame).folder, masks(iFrame + 2).name));
                    [maskCoord{iFrame,1}(:,2),maskCoord{iFrame,1}(:,1)] = ...
                        find(imgMask);
                end
            elseif smConditions.maskLoc(iMov) == 2
                imgMask = imread([qfsmPackPaths{iMov,1} fs 'QFSMPackage' fs 'refined_masks' ...
                    fs 'refined_masks_for_channel_1' fs masks(iFrame + 2).name]);
                [maskCoord{iFrame,1}(:,2),maskCoord{iFrame,1}(:,1)] = ...
                    find(imgMask);
                
            elseif smConditions.maskLoc(iMov) == 3
                imgMask = imread([qfsmPackPaths{iMov,1} fs 'SegmentationPackage' fs 'refined_masks' ...
                    fs 'refined_masks_for_channel_1' fs masks(iFrame + 2).name]);
                [maskCoord{iFrame,1}(:,2),maskCoord{iFrame,1}(:,1)] = ...
                    find(imgMask);
                
            end
        end
        
        %% Get tracks in specific FSM interval
        
        % Chop tracks into specific interval
        frameRange = [smPropPerTimeInt{iMov,1}(iFrame).smFrames(1) ...
            smPropPerTimeInt{iMov,1}(iFrame).smFrames(end)];
        
        tracksCrop = cropCompTracksFR(decompTrack,frameRange);
        
        if isempty(tracksCrop) % if this frame interval is empty (blacked out), skip to the next frame interval
            continue
        end
        
        % construct new tracksFinal from tracksCrop
        tracksFinal_new = tracksCrop;
        tracksFinal_new = rmfield(tracksFinal_new,'oldTracksInfo');
        
        %% Analyze diffusion and chop by time and span
        
        % Transient analysis & refinement (including DC-MSS) - TN20200113
        transDiffAnalysisResCrop =  basicTransientDiffusionAnalysisv1(tracksFinal_new,2,0,95); % 2nd/3rd/4th are default inputs
        
        % Chop by motion type (from DC-MSS) & by time (with MSS redone on
        % chopped tracklets)
        [refinedTrackBigSpan, measurementsBigSpan, cntTrackBigSpan] = basicTransientTrackCrop(transDiffAnalysisResCrop, tracksFinal_new); close all;
        
        % Chop by span of each tracklet within "actin domains"
        [refinedTrack, measurements, numBigSpan] = basicTransientTrackChopBySpan(refinedTrackBigSpan, measurementsBigSpan, smSpanRadius);
        
        %% Count number of tracks per chopping - Added TN20220403
        summaryNumTrack = countNumTrackletsMultChop(tracksFinal_new, transDiffAnalysisResCrop, refinedTrackBigSpan, refinedTrack, cntTrackBigSpan,numBigSpan);
        summaryNumTrack = horzcat(table(iMov, iFrame), summaryNumTrack);
        summaryNumTrackAll = vertcat(summaryNumTrackAll, summaryNumTrack);
        
        %% Get list of SM in this interval

        % reorganize cropped SM tracks information into matrix format
        trackedFeatureInfo = convStruct2MatIgnoreMS(refinedTrack,1);
        
        % matrix of x/y-coordinate of SM equipvalently ranged (col1-51)
        % appropriate to this specific frameRange.
        xyInterval = cat(3,trackedFeatureInfo(:,1 : 8 : end),trackedFeatureInfo(:,2 : 8 : end));
        
        % average position of each SM track in a particular interval
        xyMean = [nanmean(xyInterval(:,:,1),2), nanmean(xyInterval(:,:,2),2)];
        
        % enumerate SM (vector of indices of existing particles in this interval)
        indxNoNan = find(~isnan(xyMean(:,1)));
        smPropPerTimeInt{iMov,1}(iFrame,1).smList = indxNoNan;
        
        %matrix of existing SM tracks within a particular interval
        existTracks = xyInterval(indxNoNan,:,1);
                
        %% Extraction/calculation of various SM properties (many from above diffusion and assembly state analysis)
        
        % average position of SM particles relevant to this particular interval
        % stored in (x,y) coordinates
        smPropPerTimeInt{iMov,1}(iFrame,1).meanPos = [xyMean(indxNoNan,1), xyMean(indxNoNan,2)];
        
        % Amplitude of SM (matrix)
        ampMat = trackedFeatureInfo(:,4 : 8 : end);
        ampMean = nanmean(ampMat,2);
        smPropPerTimeInt{iMov,1}(iFrame,1).meanAmp = ampMean(indxNoNan);
        
        % motion classification for SM particle existing in this interval
        smPropPerTimeInt{iMov,1}(iFrame,1).trackClass = ...
            measurements(indxNoNan,1);
        
        % diffusion coefficient for SM existing in this interval
        smPropPerTimeInt{iMov,1}(iFrame,1).diffCoef = ...
            measurements(indxNoNan,3);
        
        % confinement radius
        smPropPerTimeInt{iMov,1}(iFrame,1).confRad = ...
            measurements(indxNoNan,2);
        
        % oligomeric state
        uniqueVect = nan(length(refinedTrack),1);
        for i = 1:length(refinedTrack)
            %KJ 20210715: Replaced Tra's original code (commented out)
            %with simpler line, as checking uniqueness is an overkill
            %             uniqueVect0 = unique(refinedTrack(i).aggregState);
            %             uniqueVect0 = uniqueVect0(~isnan(uniqueVect0));
            %             if length(uniqueVect0) ~= 1
            %                 error('Each decomptrack is not 1 unique oligomeric state as expected!')
            %             end
            %             uniqueVect(i,1) = uniqueVect0;
            uniqueVect(i,1) = refinedTrack(i).aggregState(1);
        end
        smPropPerTimeInt{iMov,1}(iFrame,1).aggregState = uniqueVect(indxNoNan);
        
        % calculate .meanDisp (mean displacement, for the track, over this interval),
        smPropPerTimeInt{iMov,1}(iFrame,1).meanDisp = ...
            nanmean(sqrt((xyInterval(indxNoNan,2 : end,1) - ...
            xyInterval(indxNoNan,1 : end - 1,1)) .^ 2 + ...
            (xyInterval(indxNoNan,2 : end,2) - ...
            xyInterval(indxNoNan,1 : end - 1,2)) .^ 2),2);
        
        %% SM local density
        
        availFrame = min(size(smPropPerTimeInt{iMov,1}(iFrame,1).smFrames,1),size(xyInterval,2)); % in case xyInterval doesn't have the full length of (fsmInterval+1)
        
        for iSmFrame = 1:availFrame
            
            % distance between all SM (of 1 frame)
            disMatPerFrame = createDistanceMatrix([xyInterval(:,iSmFrame,1) xyInterval(:,iSmFrame,2)], ...
                [xyInterval(:,iSmFrame,1) xyInterval(:,iSmFrame,2)]);
            
            % calculate local density of all SM (in this 1 frame)
            for iCol = 1 : size(disMatPerFrame,2) % looping through all listed SM
                
                mtchIndex = (disMatPerFrame(:,iCol) <= smConditions.rad_density(iMov));
                
                % for SM that doesn't have detection at specific frame, its
                % distance in disMatPerFrame is filled as NaN, hence its sum(mtchIndex)
                % will be 0. For all SM that exist, sum(mtchIndex) is at
                % least 1, because a detection's distance to itself is 0,
                % and it always exists inside its neighborhood.
                %
                % If SM detection exist, input density in .smDensityPerFr;
                % if SM detection doesn't exist, input NaN.
                % TN 20191010
                if  sum(mtchIndex) ~= 0
                    smPropPerTimeInt{iMov}(iFrame).smDensityPerFr(iCol,iSmFrame) = sum(mtchIndex)/...
                        (pi*((smConditions.rad_density(iMov)) ^ 2));
                else
                    smPropPerTimeInt{iMov}(iFrame).smDensityPerFr(iCol,iSmFrame) = NaN;
                end
                
            end
            
        end
        
        % calculate local SM density per FSM interval (density = number of SM per pixel)
        % .smDensity (over neighborhood of input-defined radius, averaging this interval)
        smPropPerTimeInt{iMov}(iFrame).smDensity(:,1) = ...
            nanmean(smPropPerTimeInt{iMov}(iFrame).smDensityPerFr,2);
        
        %% Calculation of various other SM properties
        
        if isempty(indxNoNan)
            
            smPropPerTimeInt{iMov,1}(iFrame,1).lifetime = zeros(0,1);
            smPropPerTimeInt{iMov,1}(iFrame,1).numObservations = zeros(0,1);
            smPropPerTimeInt{iMov,1}(iFrame,1).smFramesExist = zeros(0,1);
            smPropPerTimeInt{iMov,1}(iFrame,1).netVelocity = zeros(0,2);
            smPropPerTimeInt{iMov,1}(iFrame,1).netSpeed = zeros(0,1);
            smPropPerTimeInt{iMov,1}(iFrame,1).netVelAngle = zeros(0,1);
            smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr = zeros(0,2);
            smPropPerTimeInt{iMov,1}(iFrame,1).maskDensity = 0;
            threshMet{iMov,1}{iFrame,1} = [];
            
        else
            
            for kTrack = 1 : length(indxNoNan)
                
                % (vector of) SM frame numbers which an SM track exists in (in a particular FSM interval)
                temp = find(~isnan(existTracks(kTrack,:)))';
                
                % lifetime of this existing SM track (in a particular interval)
                smPropPerTimeInt{iMov,1}(iFrame,1).lifetime(kTrack,1) = temp(end) - temp(1) + 1;
                
                % number of observations (actual detection, not closed gap)
                % for an existing SM track (in a particular interval)
                smPropPerTimeInt{iMov,1}(iFrame,1).numObservations(kTrack,1) = ...
                    length(temp);
                
                % SM frames at which this track exists (in a particular FSM interval)
                smPropPerTimeInt{iMov,1}(iFrame,1).smFramesExist{kTrack,1} = ...
                    temp;
                                
                % Net Velocity of this existing SM track (in a particular
                % interval)
                % net velocity = (total distance displaced (|x|,|y|) in
                % pixels)/(total time passed (in number of SM frame intervals))
                smPropPerTimeInt{iMov,1}(iFrame,1).netVelocity(kTrack,:) = ...
                    (xyInterval(indxNoNan(kTrack),temp(end),:) - ...
                    xyInterval(indxNoNan(kTrack),temp(1),:)) / ...
                    (smPropPerTimeInt{iMov,1}(iFrame,1).lifetime(kTrack) - 1);
                
                % Net Speed for an existing SM track (in a particular
                % interval) (along the 2D plane)
                % net speed = sqrt(dx^2 + dy^2)/time (pixel/SM frame intervals)
                smPropPerTimeInt{iMov,1}(iFrame,1).netSpeed(kTrack,1) = ...
                    sqrt(sum(smPropPerTimeInt{iMov,1}(iFrame,1).netVelocity(kTrack,:) .^ 2));
                
                % net velocity angle (taken w.r.t image coordinate system)
                % of an existing SM track (in a particular interval)
                % range of netVelAngle = [-90, 90]
                smPropPerTimeInt{iMov,1}(iFrame,1).netVelAngle(kTrack,1) = ...
                    atand(smPropPerTimeInt{iMov,1}(iFrame,1).netVelocity(kTrack,2) / ...
                    smPropPerTimeInt{iMov,1}(iFrame,1).netVelocity(kTrack,1));
                
                % (x,y) coordinates of SM at every SM frame
                posVect = [xyInterval(indxNoNan(kTrack),:,1)' xyInterval(indxNoNan(kTrack),:,2)']; % vector of (x,y) positions of sm track
                smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack,1} = posVect(temp,:);
                
                % 20220909 TN adds randomization of SM tracks (moving SM
                % tracks from its original location to somewhere else on
                % the eroded mask)
                if isfield(smConditions,'randFlag')
                    switch smConditions.randFlag(iMov)
                        case 1 % if smConditions.randFlag(iMov) == 1, then randomize track position to any location in mask ROI
                            
                            if ismember(round(smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,:)), [maskCoord{iFrame,1}(:,1),maskCoord{iFrame,1}(:,2)],'rows')
                                % if current track is inside mask, then randomize
                                % it; if current track is outside mask, then let it
                                % be.
                                if ~isnan(smPropPerTimeInt{iMov,1}(iFrame,1).confRad(kTrack,1))
                                    % NaN span indicates that track only contain 1
                                    % detection. If that's the case, do no
                                    % randomization.
                                    
                                    % Erode current mask
                                    se_tmp = strel('disk', round(smPropPerTimeInt{iMov,1}(iFrame,1).confRad(kTrack,1)/2)); % built an element to erode by radius of SM we want to move
                                    imgMaskErode = imerode(imgMask,se_tmp);
                                    se_tmp2 = strel('disk', round(smPropPerTimeInt{iMov,1}(iFrame,1).confRad(kTrack,1)/4)); % built an element to erode by radius of SM we want to move
                                    imgMaskRefined = imopen(imgMaskErode, se_tmp2); % remove the small disconnected mask from the main mask
                                    
                                    % Find a new randomized position within eroded mask
                                    [maskCoordRefined_tmp(:,2), maskCoordRefined_tmp(:,1)] = find(imgMaskRefined);
                                    randPos_tmp = datasample(maskCoordRefined_tmp,1);
                                    %figure, imshow(imgMask), hold on, scatter(randPos_tmp(1,1),randPos_tmp(1,2),'x'); % debug to make sure randomized position is within eroded mask
                                    %scatter(smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,1),smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,2),'o'); % debug to visualize old track position (utrack & SMI-FSM report image coordinate)
                                    
                                    % Find displacement vector to translate SM track mean position to new position
                                    displacementVect_tmp = randPos_tmp - smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,:);
                                    
                                    % Translate the SM track mean position and every SM track
                                    % detections at each timepoint to new position
                                    %plot(smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack}(:,1),smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack}(:,2),'d-','Color','b') % debug to visualize old track position (image coordinate)
                                    
                                    smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,:) = randPos_tmp; % equivalent to: smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,:) + displacementVect_tmp;
                                    smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack} = smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack} + displacementVect_tmp;
                                    %plot(smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack}(:,1),smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack}(:,2),'d-','Color','y') % debug to visualize NEW track position (image coordinate)
                                    
                                    clear maskCoordRefined_tmp se_tmp se_tmp2 imgMaskErode imgMaskRefined randPos_tmp displacementVect_tmp
                                end % ( if ~isnan(smPropPerTimeInt{iMov,1}(iFrame,1).confRad(kTrack,1))  )
                            end % (   ismember(round(smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,:)), [maskCoord{iFrame,1}(:,1),maskCoord{iFrame,1}(:,2)],'rows')    )
                        
                        
                        otherwise % if smConditions.randFlag(iMov) == -M, then shift track position M pixels to the right.
                            
                            if smConditions.randFlag(iMov) < 0
                                M = abs(smConditions.randFlag(iMov)); % shift distance
                                % Translate the SM track mean position and every SM track
                                % detections at each timepoint to new position
                                smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,1) = smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,1) + M; % equivalent to: smPropPerTimeInt{iMov,1}(iFrame,1).meanPos(kTrack,:) + displacementVect_tmp;
                                smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack}(:,1) = smPropPerTimeInt{iMov,1}(iFrame,1).indiPosPerFr{kTrack}(:,1) + M;
                                
                            end % (  if smConditions.randFlag(iMov) < 0  )
                            

                    end % (   switch smConditions.randFlag(iMov)   )
                end %   (   smConditions.randFlag(iMov) == 1   )
                
            end %(for kTrack = 1 : length(indxNoNan))
            
            if smConditions.maskLoc(iMov) ~= 0
                
                % find the "valid" SM tracks whose mean position is within the mask
                threshMet{iMov,1}{iFrame,1} = find(ismember(round(smPropPerTimeInt ...
                    {iMov,1}(iFrame,1).meanPos), [maskCoord{iFrame,1}(:,1),maskCoord{iFrame,1}(:,2)],'rows'));
                
                % calculate mask density (# of tracks per pixel)
                % = number of valid SM tracks / area of mask
                smPropPerTimeInt{iMov,1}(iFrame,1).maskDensity = ...
                    length(indxNoNan(threshMet{iMov,1}{iFrame,1})) / length(maskCoord{iFrame,1}(:,1));
                
            else
                
                % if there is no mask, only find the SM that exist longer
                % than input-defined minimum lifetime
                threshMet{iMov,1}{iFrame,1} = ...
                    1:length(smPropPerTimeInt{iMov,1}(iFrame,1).lifetime);
                
                smPropPerTimeInt{iMov,1}(iFrame,1).maskDensity = NaN; %% SHOULD BE # OF TRACKS / IMAGE AREA (but there's no easy way of getting this info)
                
            end
            
        end %(if isempty(indxNoNan))
        
    end %(for iFrame = windFirst : windLast)
    
end %(for iMov = 1 : numMovs)

end

%% ~~~ the end ~~~
