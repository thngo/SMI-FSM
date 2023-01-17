function [refinedTrack, measurements, cntTrackOutput] = basicTransientTrackCrop(transDiffAnalysisRes, decompTrack)
%BASICTRANSIENTTRACKCROP splits decompTrack into tracklets/segments based on DC-MSS analysis & chop tracks in half if track duration is too long; and tally out tracks information/analysis in a matrix format
%
%SYNOPSIS: [refinedTrack, measurements] = basicTransientTrackCrop(transDiffAnalysisRes, decompTrack)
%
%INPUT:        transDiffAnalysisRes: output of DC-MSS (basicTransientDiffusionAnalysisv1.m)
%
%                       decompTrack: compound track that had DC-MSS done
%                                    on; these compound track are their own
%                                    tracklets (no merging/splitting information)
%
%OUTPUT:    refinedTrack: A structure containing "final" tracks that has only
%                         1 classification & has track size less than
%                         minTransTrackSize (<40); includes following fields:
%                                 .tracksFeatIndxCG
%                                 .tracksCoordAmpCG
%                                 .seqOfEvents
%                                 .aggregState (optional, if in input)
%
%           measurements: matrix containing information of each "final"
%                         tracks in corresponding rows of "refinedTrack":
%                         (col 1) Track Classification
%                         (col 2) Confinement Radius
%                         (col 3) Normal Diffusion Coefficient from DC-MSS
%                         "NaN"-value means that DC-MSS was not able to
%                         calculate different properties for this segment.
%
%           cntTrackOutput: (2 x 1 vector) = [inpTransTrackCnt, totTransTracklet];
%                         cntTrackOutput(1,1) reports on the number of
%                         input tracks exhibit transient motion.
%                         cntTrackOutput(1,2) reports on the number of
%                         tracklets initially come from input tracks'
%                         transient motions.
%
%Tra Ngo, Jan 2020
%
%REMARK: right now we only divide track in half (AKA max track duration =
%51 frames). Adaptations are needed in the future to make sure we devide
%into minimumDuration tracks (i.e. if we have track duration = 60 frames, it
%will be devided into 3 tracklets of 20-frame-long.)
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

%check if input contain aggregation state (which will need to be processed
%if exist)
aggregField = isfield(decompTrack,'aggregState');

%% Initialize variables
numTrack = length(transDiffAnalysisRes);
minTransTrackSize = 40; % since MSS can analyze tracks of size 20 frames, if
% a track exceeds 2x20 (minTransTrackSize), this track is chopped into 2.
mSimpleTrack = 2; % input for simpleTrackChop since track is chopped into 2.
measurements = nan(numTrack,3);

refinedTrack(numTrack,1).tracksFeatIndxCG = [];
refinedTrack(numTrack,1).tracksCoordAmpCG = [];
refinedTrack(numTrack,1).seqOfEvents = [];

if aggregField
    refinedTrack(numTrack,1).aggregState = [];
end

inpTransTrackCnt = 0; % inputTransientTrackCount (Added TN20220706)
totTransTracklet = 0; % totalTransientTracklet (Added TN20220706)

%% Loop through each (decomp)track to check if it has different transient diffusion classification by DC-MSS and process it
% STEP 1: if there are multiple classifications, chop that segment by how
% DC-MSS chops that segment.
% STEP 2: if resulted chopping has track size more than minTransTrackSize
% (>40) then chop that track in half and re-run DC-MSS or MSS.
%
% A "final" track is identified and saved in a separate "compTrack" structure
% if it has only 1 classification & has track size less than
% minTransTrackSize (<40).
% Simultaneously extract: classification, confinement radius, normal
% diffusion coefficient (MSD) from each "final" track.

for iTrack = 1:numTrack
    for iSeg = 1:length(transDiffAnalysisRes(iTrack).segmentClass)
        % If input is decompTrack, iSeg should always be 1.
        if iSeg ~= 1
            error("ERROR:: Beware that input tracks are not decompounded. Function not built for compound track in mind!");
        end
        
        %%% STEP 1: Check if number of transient classifications detected
        %%% by DC-MSS is only 1 or multiple, for this (decomp)track:
        numClass = size(transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum,1);
        
        if numClass < 1
            error("ERROR:: DC-MSS for this track is empty!")
        elseif numClass > 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Added TN20220706: To answer: How many tracks exhibit transient motion?
            inpTransTrackCnt = inpTransTrackCnt + 1; % inputTransientTrackCount
            totTransTracklet = totTransTracklet + numClass; % totalTransientTracklet
            % (to answer How many tracklets initially come from tracks exhibiting transient motion?)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % there are multiple classifications, chop segment accordingly
            % and reperform DC-MSS.
            
            currSizeRefinedTrack = size(refinedTrack,1);
            
            for iClass = 1:numClass
                absIniFrame = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(1,1);
                
                endTime = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(iClass,2);
                iniTime = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(iClass,1);
                trackClass = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(iClass,3);
                confRad = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(iClass,20);
                diffCoef = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(iClass,19);
                mssSlope = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(iClass,4); % added 20221129
                
                iniFrame = iniTime - absIniFrame + 1;
                endFrame = endTime - absIniFrame + 1;
                
                if iClass == 1
                    
                    % the first classification is copied in the same row as the original decompTrack in "refinedTrack" structure
                    
                    refinedTrack(iTrack).tracksFeatIndxCG = decompTrack(iTrack).tracksFeatIndxCG(1, iniFrame:endFrame);
                    refinedTrack(iTrack).tracksCoordAmpCG = decompTrack(iTrack).tracksCoordAmpCG(1,((iniFrame-1)*8+1):(endFrame*8));
                    
                    if aggregField
                        refinedTrack(iTrack).aggregState = decompTrack(iTrack).aggregState(1, iniFrame:endFrame);
                    end
                    
                    refinedTrack(iTrack).seqOfEvents(1,:) = decompTrack(iTrack).seqOfEvents(1,:);
                    refinedTrack(iTrack).seqOfEvents(2,1) = endTime;
                    refinedTrack(iTrack).seqOfEvents(2,2) = 2;
                    refinedTrack(iTrack).seqOfEvents(2,3) = decompTrack(iTrack).seqOfEvents(1,3);
                    refinedTrack(iTrack).seqOfEvents(2,4) = NaN;
                    
                    % Check if length exceeds pre-defined minTransTrackSize
                    
                    if (endTime - iniTime + 1) > minTransTrackSize
                        
                        %%% If sub-segment classified by DC-MSS > minTransTrackSize,
                        %%% need to chop these, run MSS.
                        %twoTrackHalved = halveTrack(refinedTrack(iTrack));
                        twoTrackHalved = simpleTrackChop(refinedTrack(iTrack),mSimpleTrack);
                        refinedTrack(iTrack) = twoTrackHalved(1);
                        refinedTrack(currSizeRefinedTrack+1) = twoTrackHalved(2);
                        
                        %%% MSS analysis on newly divided tracks
                        diffAnalysisRes1 = trackDiffusionAnalysis1(twoTrackHalved(1),1,2,0,-0.05,0,0);
                        
                        trackClass1 = diffAnalysisRes1.classification(1,2);
                        confRad1 = diffAnalysisRes1.confRadInfo.confRadius(1);
                        diffCoef1 = diffAnalysisRes1.fullDim.normDiffCoef;
                        % GET .MSSSLOPE
                        
                        diffAnalysisRes2 = trackDiffusionAnalysis1(twoTrackHalved(2),1,2,0,-0.05,0,0);
                        
                        trackClass2 = diffAnalysisRes2.classification(1,2);
                        confRad2 = diffAnalysisRes2.confRadInfo.confRadius(1);
                        diffCoef2 = diffAnalysisRes2.fullDim.normDiffCoef;
                        % GET .MSSSLOPE
                        
                        %%% Copy to measurements matrix
                        measurements(iTrack,:) = [trackClass1, confRad1, diffCoef1];
                        measurements(currSizeRefinedTrack+1,:) = [trackClass2, confRad2, diffCoef2];
                        
                    else
                        
                        %%% Copy to measurements matrix
                        measurements(iTrack,:) = [trackClass, confRad, diffCoef];
                    end
                    
                    
                else % the later classifications are copied at the bottom of "refinedTrack" structure
                    
                    currSizeRefinedTrack = size(refinedTrack,1);
                    
                    % the later classification is copied in the same row as the original decompTrack in "refinedTrack" structure
                    
                    refinedTrack(currSizeRefinedTrack+1).tracksFeatIndxCG = decompTrack(iTrack).tracksFeatIndxCG(1, iniFrame:endFrame);
                    refinedTrack(currSizeRefinedTrack+1).tracksCoordAmpCG = decompTrack(iTrack).tracksCoordAmpCG(1,((iniFrame-1)*8+1):(endFrame*8));
                    
                    if aggregField
                        refinedTrack(currSizeRefinedTrack+1).aggregState = decompTrack(iTrack).aggregState(1, iniFrame:endFrame);
                    end
                    
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(1,1) = iniTime;
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(1,2) = 1;
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(1,3) = 1;
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(1,4) = NaN;
                    
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(2,1) = endTime;
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(2,2) = 2;
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(2,3) = 1;
                    refinedTrack(currSizeRefinedTrack+1).seqOfEvents(2,4) = NaN;
                    
                    if (endTime - iniTime + 1) > minTransTrackSize
                        %%% If sub-segment classified by DC-MSS > minTransTrackSize,
                        %%% need to chop these, run MSS.
                        %twoTrackHalved = halveTrack(refinedTrack(currSizeRefinedTrack+1));
                        twoTrackHalved = simpleTrackChop(refinedTrack(currSizeRefinedTrack+1),mSimpleTrack);
                        refinedTrack(currSizeRefinedTrack+1) = twoTrackHalved(1);
                        refinedTrack(currSizeRefinedTrack+2) = twoTrackHalved(2);
                        
                        %%% MSS analysis on newly divided tracks
                        diffAnalysisRes1 = trackDiffusionAnalysis1(twoTrackHalved(1),1,2,0,-0.05,0,0);
                        
                        trackClass1 = diffAnalysisRes1.classification(1,2);
                        confRad1 = diffAnalysisRes1.confRadInfo.confRadius(1);
                        diffCoef1 = diffAnalysisRes1.fullDim.normDiffCoef;
                        
                        
                        diffAnalysisRes2 = trackDiffusionAnalysis1(twoTrackHalved(2),1,2,0,-0.05,0,0);
                        
                        trackClass2 = diffAnalysisRes2.classification(1,2);
                        confRad2 = diffAnalysisRes2.confRadInfo.confRadius(1);
                        diffCoef2 = diffAnalysisRes2.fullDim.normDiffCoef;
                        
                        %%% Copy to measurements matrix
                        measurements(currSizeRefinedTrack+1,:) = [trackClass1, confRad1, diffCoef1];
                        measurements(currSizeRefinedTrack+2,:) = [trackClass2, confRad2, diffCoef2];
                        
                    else
                        %%% Copy to measurements matrix
                        measurements(currSizeRefinedTrack+1,:) = [trackClass, confRad, diffCoef];
                        
                    end
                    
                end % ( if iClass == 1 )
            end
            
        elseif numClass == 1
            % There is only 1 classification.
            % Copy the information to "refinedTrack",
            % then check if this track has duration more than minTransTrackSize
            
            refinedTrack(iTrack) = decompTrack(iTrack);
            
            %%% STEP 2: Check if chopping has track size more than minTransTrackSize
            transTrackDuration = checkTrackDurationDCMSS(transDiffAnalysisRes(iTrack).segmentClass(iSeg));
            
            if transTrackDuration >= minTransTrackSize
                
                currSizeRefinedTrack = size(refinedTrack,1);
                
                %twoTrackHalved = halveTrack(decompTrack(iTrack));
                twoTrackHalved = simpleTrackChop(decompTrack(iTrack),mSimpleTrack);
                
                % the first halve is copied in the same row as the original decompTrack in "refinedTrack" structure
                refinedTrack(iTrack) = twoTrackHalved(1);
                % the second halve is copied at the bottom of "refinedTrack" structure
                refinedTrack(currSizeRefinedTrack+1) = twoTrackHalved(2);
                
                %%% MSS analysis on newly divided tracks
                diffAnalysisRes1 = trackDiffusionAnalysis1(twoTrackHalved(1),1,2,0,-0.05,0,0);
                
                
                trackClass1 = diffAnalysisRes1.classification(1,2);
                confRad1 = diffAnalysisRes1.confRadInfo.confRadius(1);
                diffCoef1 = diffAnalysisRes1.fullDim.normDiffCoef;
                
                diffAnalysisRes2 = trackDiffusionAnalysis1(twoTrackHalved(2),1,2,0,-0.05,0,0);
                
                trackClass2 = diffAnalysisRes2.classification(1,2);
                confRad2 = diffAnalysisRes2.confRadInfo.confRadius(1);
                diffCoef2 = diffAnalysisRes2.fullDim.normDiffCoef;
                
                
                %{
                % Previously : Re-DC-MSS analysis on newly divided tracks
                                diffAnalysisRes1T = basicTransientDiffusionAnalysisv1(twoTrackHalved(1));
                                trackClass1T = diffAnalysisRes1T.segmentClass.momentScalingSpectrum(1,3);
                                confRad1T = diffAnalysisRes1T.segmentClass.momentScalingSpectrum(1,20);
                                diffCoef1T = diffAnalysisRes1T.segmentClass.momentScalingSpectrum(1,19);
                
                                diffAnalysisRes2T = basicTransientDiffusionAnalysisv1(twoTrackHalved(2));
                                trackClass2T = diffAnalysisRes2T.segmentClass.momentScalingSpectrum(1,3);
                                confRad2T = diffAnalysisRes2T.segmentClass.momentScalingSpectrum(1,20);
                                diffCoef2T = diffAnalysisRes2T.segmentClass.momentScalingSpectrum(1,19);
                %}
                
                %%% Copy to measurements matrix
                measurements(iTrack,:) = [trackClass1, confRad1, diffCoef1];
                measurements(currSizeRefinedTrack+1,:) = [trackClass2, confRad2, diffCoef2];
                
            else
                
                %%% Copy to measurements matrix
                trackClass = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(1,3);
                confRad = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(1,20);
                diffCoef = transDiffAnalysisRes(iTrack).segmentClass(iSeg).momentScalingSpectrum(1,19);
                
                measurements(iTrack,:) = [trackClass, confRad, diffCoef];
                
            end
            
            
        end
    end % for iSeg = 1:length(transDiffAnalysisRes(iTrack).segmentClass)
end % (for iTrack = 1:numTrack)

cntTrackOutput = [inpTransTrackCnt, totTransTracklet]; % (Added TN20220706)

end



%%
function transTrackDuration = checkTrackDurationDCMSS(segmentClass)
% Check the size of the 1 segment classified by DC-MSS
% Input:        a segmentClass structure for 1 (decomp)tracklet, includes
%                   .asymmetry
%                   .momentScalingSpectrum
%                   .momentScalingSpectrum1D
endTime = segmentClass.momentScalingSpectrum(1,2);
iniTime = segmentClass.momentScalingSpectrum(1,1);
transTrackDuration = endTime - iniTime + 1;
end


% % % %% Obsolete sub-function (being replaced by simpleTrackChop.m)
% % % function twoTrackHalved = halveTrack(decompTrack)
% % % % Chop the decompTrack into half and return the two halves.
% % % % Input: decompTrack should only contain 1 row (1 track)
% % % % decompTrack and twoTrackHalved are structure with fields
% % % %   .tracksFeatIndxCG
% % % %   .tracksCoordAmpCG
% % % %   .seqOfEvents
% % % twoTrackHalved(2,1).tracksFeatIndxCG = [];
% % % twoTrackHalved(2,1).tracksCoordAmpCG = [];
% % % twoTrackHalved(2,1).seqOfEvents = [];
% % %
% % % % size of chopping
% % % divTransTrackSize = floor(size(decompTrack.tracksFeatIndxCG,2)/2);
% % %
% % % twoTrackHalved(1).tracksFeatIndxCG = decompTrack.tracksFeatIndxCG(1, 1:divTransTrackSize);
% % % twoTrackHalved(1).tracksCoordAmpCG = decompTrack.tracksCoordAmpCG(1,(1:(divTransTrackSize*8)));
% % %
% % % twoTrackHalved(2).tracksFeatIndxCG = decompTrack.tracksFeatIndxCG(1, (divTransTrackSize+1):end);
% % % twoTrackHalved(2).tracksCoordAmpCG = decompTrack.tracksCoordAmpCG(1,((divTransTrackSize)*8+1):end);
% % %
% % % % create new .seqOfEvents
% % % twoTrackHalved(1).seqOfEvents(1,:) = decompTrack.seqOfEvents(1,:);
% % %
% % % twoTrackHalved(1).seqOfEvents(2,1) = decompTrack.seqOfEvents(1,1) + divTransTrackSize - 1;
% % % twoTrackHalved(1).seqOfEvents(2,2) = 2;
% % % twoTrackHalved(1).seqOfEvents(2,3) = decompTrack.seqOfEvents(1,3);
% % % twoTrackHalved(1).seqOfEvents(2,4) = NaN;
% % %
% % % twoTrackHalved(2).seqOfEvents(1,1) = twoTrackHalved(1).seqOfEvents(2,1) + 1;
% % % twoTrackHalved(2).seqOfEvents(1,2) = 1;
% % % twoTrackHalved(2).seqOfEvents(1,3) = decompTrack.seqOfEvents(1,3);
% % % twoTrackHalved(2).seqOfEvents(1,4) = NaN;
% % %
% % % twoTrackHalved(2).seqOfEvents(2,:) = decompTrack.seqOfEvents(2,:);
% % % end