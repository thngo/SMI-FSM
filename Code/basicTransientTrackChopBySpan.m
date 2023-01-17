function [refinedTrackSpan, measurementsSpan, numBigSpan] = basicTransientTrackChopBySpan(refinedTrack, measurements, spanRadius)
%BASICTRANSIENTTRACKCHOPBYSPAN divide tracklets (from refinedTrack) that span long distance into spatially localized tracklets within an input acceptable span (by max pairwise distance) and tally out tracks information/analysis in a matrix format
%
%SYNOPSIS: [refinedTrackSpan, measurementsSpan] = basicTransientTrackChopBySpan(refinedTrack, measurements, spanRadius)
%
%INPUT:        refinedTrack: output of frame-based tracklets refinement
%                            (basicTransientTrackCrop.m)
%                            (structure) contains "final" tracks that has only
%                            1 classification (homogeneous diffusion) & has track
%                            size less than minTransTrackSize (<40);
%                            includes following fields:
%                                 .tracksFeatIndxCG
%                                 .tracksCoordAmpCG
%                                 .seqOfEvents
%                                 .aggregState (optional, if in input)
%
%              measurements: output of frame-based tracklets refinement
%                            (basicTransientTrackCrop.m)
%                           (matrix) contains information of each "final"
%                         tracks in corresponding rows of "refinedTrack":
%                         (col 1) Track Classification
%                         (col 2) Confinement Radius
%                         (col 3) Normal Diffusion Coefficient from DC-MSS
%
%                spanRadius: (number) radius such that if tracklets' radius
%                           is ~M*spanRadius, tracklets would be splited
%                           into M tracklets.
%                           Optional. Default: 11.1 (pixel, radius of spatially localized domain)
%
%OUTPUT:    refinedTrackSpan: (structure) contains "final" tracks that has only
%                         1 classification spanning within input spanRadius
%                         includes following fields:
%                                 .tracksFeatIndxCG
%                                 .tracksCoordAmpCG
%                                 .seqOfEvents
%
%           measurementsSpan: (matrix) contains information of each "final"
%                         tracks in corresponding rows of "refinedTrack_span":
%                         (col 1) Track Classification
%                         (col 2) Confinement Radius
%                         (col 3) Normal Diffusion Coefficient from DC-MSS
%                         "NaN"-value means that DC-MSS was not able to
%                         calculate different properties for this segment.
%
%           numBigSpan      : (integer) number of input tracks that has big
%                             span and needs to be chopped
%
%           REMARK: information and measurement of chopped-out tracklets is
%           appended to the end of the refinedTrackSpan structure and the
%           bottom of the measurementsSpan matrix.
%
%Tra Ngo, Mar 2020
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

%% Control input

if length(refinedTrack) ~= size(measurements,1)
    error('refinedTrack and measurements must have matching number of tracks (rows)');
end

if nargin < 3
    spanRadius = 300/9; % 3 micrometer is an acceptable domain to average speckle properties.
end

%% Initialize output variable & const

colDiffRad = 2; % diffusion radii of tracks are stored in col#2 of input "measurements"
refinedTrackSpan = refinedTrack;
measurementsSpan = measurements;
% Input parameters for estimConfRad
%probDim = 2; % default from basicTransientDiffusionAnalysisv1.m (problem dimension)
%confRadMin = 0; % default from estimConfRad.m (mean(eigenVal) is used)

%% Calculate span for all input tracks as MPD in each SM (remark: measurementsSpan now differs from measurements)

diffRadOut = getTrackMpd(refinedTrack); % [diffRadOut, detectPairOfMpd, mpdOriginal] = getTrackMpd(refinedTrack);
% replace the confinement radius with span (MPD)
measurementsSpan(:,colDiffRad) = diffRadOut;

%% Find the tracklets that span long distance > spanRadius

indxBigSpan = find(measurementsSpan(:,colDiffRad) > spanRadius);
numBigSpan = size(indxBigSpan,1);

%% Splitting into M tracklets
while any(indxBigSpan)
    
    % Number of m tracklets a single simple tracklet is to be split into
    numSplit = ceil(measurementsSpan(indxBigSpan,colDiffRad)./spanRadius);
    
    % REMARK: Since each track has homogeneous diffusion properties, if a tracklet
    % spans > m*spanRadius (spatially localized domain), then we split tracks
    % into m sub-tracklets. (i.e. chop into 2 if spanning 2 domains, into 3
    % if spanning 3 domains, etc.)
    
    for iIndx = 1:length(indxBigSpan)
        indxCurrent = indxBigSpan(iIndx);
        numMax = length(refinedTrackSpan);
        
        mTrackChopped = simpleTrackChop(refinedTrackSpan(indxCurrent), numSplit(iIndx)); % structure ARRAY (not structure, each entry is a structure)
                
        % First row of mTrackChopped replaces the original track in refinedTrack
        refinedTrackSpan(indxCurrent) = mTrackChopped(1);
        
        % update new span (MPD)
        spanTmp = getTrackMpd(refinedTrackSpan(indxCurrent));
        measurementsSpan(indxCurrent,colDiffRad) = spanTmp;
        
        % Subsequent rows of mTrackChopped are added to the end of refinedTrack
        for iAdd = 1:(max(size(mTrackChopped))-1) % length(mTrackChopped)-1 b/c length will give the biggest dimension % sometimes structure is saved as 1 x 2, sometimes 2 x 1
            
            refinedTrackSpan(numMax+iAdd) = mTrackChopped(iAdd+1);
            measurementsSpan(numMax+iAdd,:) = measurementsSpan(indxCurrent,:); % copy in new row (of new chopped tracklet) properties from old "mother" track
            
            % update new span (MPD)
            spanTmp = getTrackMpd(refinedTrackSpan(numMax+iAdd));
            measurementsSpan(numMax+iAdd, colDiffRad) = spanTmp;
            
        end
        
    end % (iIndx = 1:length(indxBigSpan))
    
    indxBigSpan = find(measurementsSpan(:,colDiffRad) > spanRadius);
    
end % (while any(indxBigSpan))

end