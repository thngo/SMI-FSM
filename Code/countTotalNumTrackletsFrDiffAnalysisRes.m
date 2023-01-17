function [numTrack, numAdditionalTrack] = countTotalNumTrackletsFrDiffAnalysisRes(transDiffAnalysisResCrop)
%COUNTTOTALNUMTRACKLETS look at output of DC-MSS and count how many tracklets there are in total
%
%INPUT: transDiffAnalysisResCrop (cell array) output of
%basicTransientDiffusionAnalysisv1.m
%
%OUTPUT:           numTrack: (integer) number of tracklets categorized by DC-MSS
%
%        numAdditionalTrack: (integer) number of additional tracklets if we
%                             separate the tracklets by DC-MSS
%                          numAdditionalTrack = numTrack - numOriginalTrack;
%
%Tra Ngo, Apr 2022
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
numOriginalTrack = size(transDiffAnalysisResCrop,1);
numTrackVect = nan(numOriginalTrack ,1);
for iTrack = 1:length(transDiffAnalysisResCrop)
    for iSeg = 1:length(transDiffAnalysisResCrop(iTrack).segmentClass)
        % If input is decompTrack, iSeg should always be 1.
        if iSeg ~= 1
            error("ERROR:: Beware that input tracks are not decompounded. Function not built for compound track in mind!");
        end
        % concaternate number of tracklet from this original full track
        numClass = size(transDiffAnalysisResCrop(iTrack).segmentClass(iSeg).momentScalingSpectrum,1);
        numTrackVect(iTrack) = numClass;
    end
end

numTrack = sum(numTrackVect);
numAdditionalTrack = numTrack - numOriginalTrack;