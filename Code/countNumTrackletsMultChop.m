function summaryNumTrackTab  = countNumTrackletsMultChop(tracksFinal, transDiffAnalysis, refinedTrackBigSpan, refinedTrackSmallSpan, cntTrackBigSpan, numBigSpan)
%COUNTNUMTRACKLETSMULTCHOP summarize number of tracklets added by each chopping
%
%INPUT:  tracksFinal          : output of u-track, decompounded
%        transDiffAnalysis    : output of DC-MSS
%        refinedTrackBigSpan  : output of chopping by DC-MSS and chopping
%                               in time, decompounded 
%        refinedTrackSmallSpan: output of chopping by span
%        cntTrackBigSpan      : (1x2 vector), output of chopping by DC-MSS with
%                               function basicTransientTrackCrop.m
%
%        Example input: 
%        countNumTrackletsMultChop(tracksFinal_new, transDiffAnalysisResCrop, refinedTrackBigSpan, refinedTrack, cntTrackBigSpan)
%
%OUTPUT: summaryNumTrackTab: (table) number of tracklets added by each chopping
%                           
% Tra Ngo, Apr, 2022
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

%% Number of tracklets (decompounded tracks) before chopping
numTrack = size(tracksFinal,1);
% disp(['Number of tracklets before any chopping: ' num2str(numTrack) '.'])

%% Number of tracklets 
[numTrackDcmss, addedByDcmss] =  countTotalNumTrackletsFrDiffAnalysisRes(transDiffAnalysis);
%disp(['Number of tracklets after DC-MSS: ' num2str(numTrackDcmss) '.'])
%disp(['Number of tracklets added by chopping by DC-MSS: ' num2str(addedByDcmss) '.'])

%% Number of tracklets after chopping by time:
numAfterTime = size(refinedTrackBigSpan,1);
%disp(['Number of tracklets after chopping by time: ' num2str(numAfterTime) '.'])
addedByTime = numAfterTime - numTrackDcmss;
%disp(['Number of tracklets added by chopping by time: ' num2str(addedByTime) '.'])

%% Number of tracklets after chopping by span:
numAfterSpan = size(refinedTrackSmallSpan,1);
addedBySpan = numAfterSpan - numAfterTime;
%disp(['Number of tracklets after chopping by time: ' num2str(numAfterSpan) '.'])
%disp(['Number of tracklets added by chopping by time: ' num2str(addedBySpan) '.'])

trackWithTransMot = cntTrackBigSpan(1,1);
totTransMotTracklet = cntTrackBigSpan(1,2);

summaryNumTrackTab = table(numTrack,trackWithTransMot, totTransMotTracklet, ...
    addedByDcmss,numTrackDcmss,addedByTime,numAfterTime,addedBySpan,numAfterSpan, ...
    numBigSpan);
end
