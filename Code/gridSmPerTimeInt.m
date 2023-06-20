function [smPropPerTimeBoxed, smTrackletIdPerBox] = gridSmPerTimeInt(smPropPerTimeInt, varargin)
%GRIDSMPROPPERTIMEINT: split movie frame into non-overlapping boxes and get SM information within each box
%
%SYNOPSIS: [smPropPerTimeBoxed, smTrackletIdPerBox] = gridSmPerTimeInt(smPropPerTimeInt, imageSize, varargin)
%
%EXAMPLE: smPropPerTimeBoxed = gridSmPropPerTimeInt(smPropPerTimeIntSpan100,'imageSize', [341 800]);
%
%INPUT:
%Required:
%           smPropPerTimeInt: (n x 1 cell array of structures) SMI arm
%                             analysis results of SMI-FSM analysis for n
%                             movies; each structure contain SM properties
%                             for m FSM intervals. For details,
%                             see smPropPerTimeInterval.m
%
%Optional: (as name-value pair)
%           boxWidth  : box width (AKA. number of pixel for each box
%                       row-wise)
%                       Default: 15
%
%           boxHeight : box height (AKA. number of pixel for each box
%                       column-wise)
%                       Default: 15
%
%           nBorder   : number of pixel away from image border that we are
%                       considering the image. This is to correct for dark
%                       pixel around the image border.
%                       Default: 0
%
%           imageSize : image size, determined by the microscope system
%                       (i.e., [341 800] for Olympic Scope, [512 512] for
%                       Nikon Scope, programed on Metomorph).
%                       This can be determined by MovieData.imSize_
%                       Default: [512 512]
%
%           maskInfo  : (array of cells) each cell is either a path to
%                       where the mask is stored, or is the mask matrix
%                       itself.
%                       For folder provided, load each mask per frame.
%                       For matrix provided, assume 1 mask for all frames of 1 cell.
%                       Default: [] (no mask, take everything boxes in grid)
%
%OUTPUT:   smPropPerTimeBoxed: (cell array of structure) similar to the
%                               input (smPropPerTimeInt) but with an
%                               extra field '.gridPosition' indicating the
%                               box index that the SM fall into within
%                               the grid.
%
%          smTrackletIdPerBox : (cell array of structure) ID of SM per box,
%                          number of box is dependent on imageSize and grid
%                          size. Each cell contains the following field:
%                       .boxId         : list the linear index of the boxes,
%                                        from 1 to N, N is dependent on
%                                        image size and box size.
%                       .trackletId    : number ID of tracklets that falls
%                                        within a particular box. ID is
%                                        from listing in smPropPerTimeInt.
%                       .timeRelToStart: time relative to start. This
%                                        is tracklet localizations in order
%                                        of first to last localization.
%                       .boxSize       : the number of pixel within a
%                                        particular box. This is dependent
%                                        on image size and box size. Boxes
%                                        at the image edge will be smaller
%                                        if image size is not divisible by
%                                        box size.
%
% Huong Tra Ngo, Jan 2023
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
% along with SMI-FSM. If not, see <http://www.gnu.org/licenses/&gt;.
% 

%% Input
ip = inputParser;
ip.CaseSensitive = false;
ip.FunctionName = 'gridSmPropPerTimeInt';

% Required
ip.addRequired('smPropPerTimeInt', @iscell) % expect actinPropPerTimeInt to be cell array, each cell correspond to 1 movie

% Optional
ip.addParameter('boxWidth', 15, @isscalar)
ip.addParameter('boxHeight', 15, @isscalar)
ip.addParameter('nBorder', 0, @isscalar)
ip.addParameter('imageSize', [512 512], @isvector)
ip.addParameter('maskInfo', [], @iscell) % mask location or mask; if box is outside of mask => value is NaN

ip.parse(smPropPerTimeInt, varargin{:})

% (Convenience)
boxWidth = ip.Results.boxWidth;
boxHeight = ip.Results.boxHeight;
nBorder = ip.Results.nBorder;
imageSize = ip.Results.imageSize;
maskInfo = ip.Results.maskInfo;

numMov = length(smPropPerTimeInt);

mockBox = boxImage(nan(imageSize), [boxWidth boxHeight], 'nanmean'); % create a tmp boxed image
lenBoxList = length(mockBox(:));
clear mockBox

se = strel("rectangle", [boxWidth boxHeight]);
numErosion = 3; % number of iterative erosion done on mask

%% Output
smPropPerTimeBoxed = smPropPerTimeInt;
smTrackletIdPerBox = cell(numMov,0);

%% Mask
numFr = nan(numMov,1);
imgMask = cell(numMov,0);

for kMov = 1:numMov
    numFr(kMov) = size(smPropPerTimeInt{kMov,1},1);
    for kFr = 1:numFr(kMov)
        
        if ischar(maskInfo{kMov})
            % If mask locations are provided, then go to that folder and load the masks.
            % i.e. maskInfo{1} = /project/biophysics/jaqaman_lab/vegf_tsp1/adasgu/FSM_SMI_ABD_ABDRA/20200113_PD80_LP70_Transfex_point25ug_CD36_C2_D2_day2_halo_3nMJF549_TIMEmngActlowP14/CD36/Actin/m-01-02/QFSMPackage/refined_masks/refined_masks_for_channel_1
            tmp = dir(maskInfo{kMov});
            tmp2 = imread([maskInfo{kMov} filesep tmp(kFr + 2).name]);
            imgMask{kMov,1}{kFr,1} = double(tmp2);
            imgMask{kMov,1}{kFr,1}(imgMask{kMov,1}{kFr,1}==0) = NaN;
            
        else
            disp('Program assumes user inputs masks for maskInfo');
        end
        
        % Sense check: image mask should have the same size as input mask image
        if any(size(imgMask{kMov,1}{kFr,1}) ~= imageSize)
            error('Loaded image has a different size than input image size')
        end
        
    end % for kFr = 1: length(actinPropPerTimeInt{kMov,1})
end % parfor kMov = 1:numMov

clear kMov kFr

%% Grid each movie
for i = 1:numMov
    if isfield(smPropPerTimeInt{i,1},'mergeInfoSpace') && isfield(smPropPerTimeInt{i,1},'splitInfoSpace') && isfield(smPropPerTimeInt{i,1},'smFrames') && isfield(smPropPerTimeInt{i,1},'maskDensity')
        smStruct = rmfield(smPropPerTimeInt{i,1},{'smFrames','maskDensity','mergeInfoSpace', 'splitInfoSpace'}); % remove unnecessary fields that do not describe individual tracklets
    else
        smStruct = smPropPerTimeInt{i,1};
    end
    numFr = size(smStruct,1);
    
    smTrackletIdPerBox{i,1} = repmat(struct('boxId',[],'trackletId',{},'timeRelToStart',{}), numFr,1);
    
    for iFr = 1:numFr  % looping through each FSM frame
        trackNum = size(smStruct(iFr).indiPosPerFr,1);
        % Initialize new field:
        smPropPerTimeBoxed{i,1}(iFr).gridPosition = cell(trackNum,1);
        
        smTrackletIdPerBox{i,1}(iFr,1).boxId = [1:lenBoxList]';
        smTrackletIdPerBox{i,1}(iFr,1).trackletId = cell(lenBoxList,1);
        smTrackletIdPerBox{i,1}(iFr,1).timeRelToStart = cell(lenBoxList,1);
        smTrackletIdPerBox{i,1}(iFr,1).mergePosition = cell(lenBoxList,1);
        smTrackletIdPerBox{i,1}(iFr,1).splitPosition = cell(lenBoxList,1);
        
        % STEP 0: If mask is provided, process it similarly (crop corner + box)
        % and incorporate mask information. (NaN out boxes that are outside of masks)
        
        if ~isempty(imgMask)
            % cropImgMask = imgMask{i,1}{iFr,1}(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
            
            % Create map that categorize boxes that are within/partially within/outside mask:
            numErosion = 0; se = strel('rectangle', [boxWidth boxHeight]); senseCheck = false;
            [~, gridClassInfo] = classifyMaskGridByErosion(imgMask{i,1}{iFr,1}, [boxWidth boxHeight], numErosion, se, nBorder, senseCheck);
            
            for iBox = 1:lenBoxList
                smTrackletIdPerBox{i,1}(iFr,1).boxSize{iBox,1} = gridClassInfo.boxSize(iBox);
                smTrackletIdPerBox{i,1}(iFr,1).maskSize{iBox,1} = gridClassInfo.maskSize(iBox);
                smTrackletIdPerBox{i,1}(iFr,1).maskStatus{iBox,1} = gridClassInfo.maskStatus(iBox);
            end % ( for iBox = 1:lenBoxList )
            
        end % ( if ~isempty(imgMask) )
        
        for iTrk = 1:trackNum % Looping through each valid SM tracks
            % STEP 1: Construct image with enumeration at each SM localization:
            [constructImg, constructImgArray] = convCoord2Img(round(smStruct(iFr).indiPosPerFr{iTrk, 1}), imageSize, 'enumerate', 'nan');
            
            outputImg = constructImg(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
            
            % STEP 2: Grid the reconstructed image:
            boxList = boxImage(outputImg, [boxWidth boxHeight], 'nanmean'); % create a tmp boxed image
            
            % STEP 3: index back into smPropPerTimeInterval the box that SM
            % localization of the current track sit in:
            for iBox = 1:lenBoxList
                if any(any(~isnan(boxList{iBox}))) % Boxes outside the mask are indicated as NaN and are not considered.
                    tmp = ~isnan(boxList{iBox});
                    indexTmp = boxList{iBox}(tmp);
                    smPropPerTimeBoxed{i,1}(iFr,1).gridPosition{iTrk, 1}(indexTmp,1) = iBox;
                    
                    % STEP 4: Update boxes intentified here with corresponding SM identity:
                    smTrackletIdPerBox{i,1}(iFr,1).trackletId{iBox,1} = vertcat(smTrackletIdPerBox{i,1}(iFr,1).trackletId{iBox,1}, iTrk);
                    smTrackletIdPerBox{i,1}(iFr,1).timeRelToStart{iBox,1} = [smTrackletIdPerBox{i,1}(iFr,1).timeRelToStart{iBox,1}; {indexTmp}];
                    smTrackletIdPerBox{i,1}(iFr,1).boxSize{iBox,1} = length(boxList{iBox}(:));
                    
                    % Collect the boxes with any SM tracks in it
                    % for j = 1:length(smTrackletIdPerBox{1,1}(iFr).trackId), if ~isempty(smTrackletIdPerBox{1,1}(iFr).trackId{j}), disp(['!!!!' num2str(j)]); end; end
                    
                end
            end % ( for iBox = 1:lenBoxList )
            
            
            for iImage = 1:length(constructImgArray) % if constructImgArray is empty, meaning there is no pixel with repeated localization, then this for-loop is never entered.
                
                outputImg = constructImgArray{iImage}(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
                
                % STEP 2: Grid the reconstructed image:
                boxList = boxImage(outputImg, [boxWidth boxHeight], 'nanmean'); % create a tmp boxed image
                
                % STEP 3: index back into smPropPerTimeInterval the box that SM
                % localization of the current track sit in:
                for iBox = 1:lenBoxList
                    if any(any(~isnan(boxList{iBox}))) % Boxes outside the mask are indicated as NaN and are not considered.
                        tmp = ~isnan(boxList{iBox});
                        indexTmp = boxList{iBox}(tmp); % NEED TO MAKE THIS NAN BECAUSE (1) THE REST OF INDEXTMP THAT DOES NOT APPEAR BECOME ASSIGNED TO BOX 0
                        % WHICH CAN BE CONFUSING AND
                        % WHEN I CHECK THE SM TRACK ON CELL MASK, GRID POSITION DOES NOT SEEM CORRECTLY ASSIGNED
                        smPropPerTimeBoxed{i,1}(iFr,1).gridPosition{iTrk, 1}(indexTmp,1) = iBox;
                        
                        % STEP 4: Update boxes intentified here with corresponding SM identity:
                        smTrackletIdPerBox{i,1}(iFr,1).trackletId{iBox,1} = vertcat(smTrackletIdPerBox{i,1}(iFr,1).trackletId{iBox,1}, iTrk);
                        smTrackletIdPerBox{i,1}(iFr,1).timeRelToStart{iBox,1} = [smTrackletIdPerBox{i,1}(iFr,1).timeRelToStart{iBox,1}; {indexTmp}];
                        smTrackletIdPerBox{i,1}(iFr,1).boxSize{iBox,1} = length(boxList{iBox}(:));
                        
                        % 1-liner to collect the boxes with any SM tracks in it
                        % for j = 1:length(smTrackletIdPerBox{1,1}(iFr).trackId), if ~isempty(smTrackletIdPerBox{1,1}(iFr).trackId{j}), disp(['!!!!' num2str(j)]); end; end
                        
                    end
                end % ( for iBox = 1:lenBoxList )
                
            end % for iImage = 1:length(constructImgArray)
            
            % STEP 5: Get the 1 box that track this matched to most of the time
            uniqueBox = unique(smPropPerTimeBoxed{i,1}(iFr,1).gridPosition{iTrk, 1});
            numLocalization = length(smPropPerTimeBoxed{i,1}(iFr,1).gridPosition{iTrk, 1});
            percInBox = nan(length(uniqueBox),1);
            for iUniq = 1: length(uniqueBox)
                percInBox(iUniq) = nansum(smPropPerTimeBoxed{i,1}(iFr,1).gridPosition{iTrk, 1} == uniqueBox(iUniq))/numLocalization;
            end
            indxBoxAssigned = percInBox == max(percInBox);
            smPropPerTimeBoxed{i,1}(iFr,1).boxAssigned{iTrk, 1} = uniqueBox(indxBoxAssigned);
            smPropPerTimeBoxed{i,1}(iFr,1).percInBoxAssigned{iTrk, 1} = percInBox(indxBoxAssigned);
            
        end % ( for iTrk = 1:trackNum )
        
        %% Get splitting & merging information - TN20230324:
        if isfield(smPropPerTimeInt{i,1},'mergeInfoSpace') && isfield(smPropPerTimeInt{i,1},'splitInfoSpace')
            mergeLoc = smPropPerTimeInt{i,1}.mergeInfoSpace;
            splitLoc = smPropPerTimeInt{i,1}.splitInfoSpace;
            
            % STEP 1: Construct image with enumeration at each MS localization:
            mergeImg = convCoord2Img(round(mergeLoc), imageSize, 'enumerate', 'nan');
            mergeImgOut = mergeImg(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
            splitImg = convCoord2Img(round(splitLoc), imageSize, 'enumerate', 'nan');
            splitImgOut = splitImg(1 + nBorder:end - nBorder, 1 + nBorder:end - nBorder);
            
            % STEP 2: Grid the reconstructed image:
            boxListM = boxImage(mergeImgOut, [boxWidth boxHeight], 'nanmean'); % create a tmp boxed image
            boxListS = boxImage(splitImgOut, [boxWidth boxHeight], 'nanmean'); % create a tmp boxed image
            
            % STEP 3: index back into smPropPerTimeInterval the box that SM
            % localization of the current track sit in:
            for iBox = 1:lenBoxList
                if any(any(~isnan(boxListM{iBox})))
                    tmp = ~isnan(boxListM{iBox});
                    indexTmp = boxListM{iBox}(tmp);
                    smTrackletIdPerBox{i,1}(iFr,1).mergePosition{iBox,1} = vertcat(smTrackletIdPerBox{i,1}(iFr,1).mergePosition{iBox,1}, mergeLoc(indexTmp,:));
                end
                
                if any(any(~isnan(boxListS{iBox})))
                    tmp = ~isnan(boxListS{iBox});
                    indexTmp = boxListS{iBox}(tmp);
                    smTrackletIdPerBox{i,1}(iFr,1).splitPosition{iBox,1} = vertcat(smTrackletIdPerBox{i,1}(iFr,1).splitPosition{iBox,1}, splitLoc(indexTmp,:));
                end
            end % ( for iBox = 1:lenBoxList )
            
        end
    end % ( for iFr = 1: size(smPropStruct,1)  )
end % ( parfor i = 1:numMov )

end