function [decompTrack] = decompoundCompTracks(compTrack)
%DECOMPOUNDCOMPTRACKS decompound (separate) all interacting segments of comptracks into separate tracks (ignore merge-split/MS)
%
%SYNOPSIS:  [decompTrack] = decompoundCompTracks(compTrack)
%
%INPUT   :  compTrack : compound tracks structure (in alternative format),
%                       each track can include multiple interacting
%                       segments; includes following fields:
%                  .tracksFeatIndxCG
%                  .tracksCoordAmpCG
%                  .seqOfEvents: some track segments (tracklets) may not
%                  have initiating/terminating time point as input was
%                  processed through alterative format.
%                  , etc.
%                  .aggregState: (optional) same dimensions as tracksFeatIndxCG
%
%OUTPUT  :  decompTrack: compound tracks structure, each track only
%                        contains 1 segment; includes following fields:
%                  .tracksFeatIndxCG
%                  .tracksCoordAmpCG
%                  .seqOfEvents: all tracklet starts and ends become "true" events
%                  (AKA we lose all merging & splitting event information)
%                  .aggregState (optional, only if input contains this
%                  field)
%
%Remark: does not check for starting or running (trailing) 0s. (Mar 2021)
%Tra Ngo, Jan 2020. Add handling of .aggregState in Sep 2020
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

%% Check input
if isempty(compTrack) || ~isfield(compTrack,"tracksFeatIndxCG") || ~isfield(compTrack,"tracksCoordAmpCG") || ~isfield(compTrack,"seqOfEvents")
    error("Input compound track format invalid!");
end

%check if input contain aggregation state (which will need to be processed
%if exist)
aggregField = isfield(compTrack,'aggregState');

%% Initiate new variable
numCompTrack = length(compTrack);

% Initiate decompTrack structure with at least same number of entries as compTrack
decompTrack(numCompTrack,1).tracksFeatIndxCG = [];
decompTrack(numCompTrack,1).tracksCoordAmpCG = [];
decompTrack(numCompTrack,1).seqOfEvents = [];

if aggregField
    decompTrack(numCompTrack,1).aggregState = [];
end

iDecompTrack = 1;

%% Fill in decompTrack structure by going through each track of compTrack
for iCompTrack = 1:numCompTrack % loop through each compTrack
    
    numTracklet = size(compTrack(iCompTrack).tracksFeatIndxCG,1);
    
    if numTracklet == 1 %if compound track is really only 1 tracklet, simply copy it
        
        decompTrack(iDecompTrack).tracksFeatIndxCG = compTrack(iCompTrack).tracksFeatIndxCG;
        decompTrack(iDecompTrack).tracksCoordAmpCG = compTrack(iCompTrack).tracksCoordAmpCG;
        decompTrack(iDecompTrack).seqOfEvents = compTrack(iCompTrack).seqOfEvents;
        if aggregField
            decompTrack(iDecompTrack).aggregState = compTrack(iCompTrack).aggregState;
        end
        
        %%% increment iDecompTrack after filling in each row of decompTrack structure
        iDecompTrack = iDecompTrack+1;
        
    else %if compound track contains more than 1 tracklet
        
        %%% Determine time frame of .seqOfEvents
        iniFrame = compTrack(iCompTrack).seqOfEvents(1,1);
        endFrame = compTrack(iCompTrack).seqOfEvents(end,1);
        totalInt = endFrame - iniFrame + 1;
        
        for iTracklet = 1:numTracklet % loop through each tracklet inside compTrack
            
            indxTracklet = compTrack(iCompTrack).seqOfEvents(:,3) == iTracklet;
            
            % Copy seqOfEvents of only tracklet of interest
            decompTrack(iDecompTrack).seqOfEvents = compTrack(iCompTrack).seqOfEvents(indxTracklet,:);
            
            % some track segments (tracklets) may not have initiating/terminating
            % (or neither) time points as input was processed through alterative format.
            % If so, create .seqOfEvents for them
            
            % Note: if both initiating and terminating time points exist, no
            % case in the switch-case is entered.
            
            switch sum(indxTracklet)
                
                case 1 % only either initiating or terminating time point exist
                    
                    if decompTrack(iDecompTrack).seqOfEvents(1,2) == 1 % only Initiating time exists => tracklet terminates as
                        % other tracklet merges with it => search for terminating time
                        
                        indxTrackletTerm = ... % find rows where tracklet could potentially terminate
                            find(compTrack(iCompTrack).seqOfEvents(:,4) == decompTrack(iDecompTrack).seqOfEvents(1,3));
                        
                        if length(unique(compTrack(iCompTrack).seqOfEvents(indxTrackletTerm,1))) == 1
                            decompTrack(iDecompTrack).seqOfEvents(2,:) = compTrack(iCompTrack).seqOfEvents(indxTrackletTerm(1),:);
                        else
                            error("ERROR: Tracklet terminates at multiple time points!");
                        end
                        
                    elseif decompTrack(iDecompTrack).seqOfEvents(1,2) == 2 % only Terminating time exists => search for initiating time
                        
                        indxTrackletIni = ...
                            find(compTrack(iCompTrack).seqOfEvents(:,4) == decompTrack(iDecompTrack).seqOfEvents(1,3));
                        decompTrack(iDecompTrack).seqOfEvents(2,:) = decompTrack(iDecompTrack).seqOfEvents(1,:); % shift 1st row to 2nd row
                        
                        if length(unique(compTrack(iCompTrack).seqOfEvents(indxTrackletIni,1))) == 1 % only 1 time point detected
                            decompTrack(iDecompTrack).seqOfEvents(1,:) = compTrack(iCompTrack).seqOfEvents(indxTrackletIni(1),:);
                        else
                            error("ERROR: Tracklet initiates at multiple time points!");
                        end
                        
                    end
                    
                case 0 % neither initiating nor terminating time point exist
                    
                    % search for initiating time
                    indxTrackletIni = ...
                        find(compTrack(iCompTrack).seqOfEvents(:,4) == iTracklet & ...
                        compTrack(iCompTrack).seqOfEvents(:,2) == 2);
                    
                    % search for terminating time
                    indxTrackletTerm = ...
                        find(compTrack(iCompTrack).seqOfEvents(:,4) == iTracklet & ...
                        compTrack(iCompTrack).seqOfEvents(:,2) == 1);
                    
                    % copy .seqOfEvents if initiating time < terminating time
                    if compTrack(iCompTrack).seqOfEvents(indxTrackletIni(1),1) < compTrack(iCompTrack).seqOfEvents(indxTrackletTerm(1),1)
                        decompTrack(iDecompTrack).seqOfEvents(1,:) = compTrack(iCompTrack).seqOfEvents(indxTrackletIni(1),:);
                        decompTrack(iDecompTrack).seqOfEvents(2,:) = compTrack(iCompTrack).seqOfEvents(indxTrackletTerm(1),:);
                    else
                        error("ERROR: Initiating time > Terminating time!")
                    end
                    
            end %(switch sum(indxTracklet))
            
            %% *Re-format Notes:: seqOfEvent so that all events are true (col4 = NaN)
            % and segment number is similar between row indicating start time
            % and end time.
            % Note: When a tracklet ends by transition into another tracklet,
            % no detection in observed for this tracklet at this "terminating" frame.
            % So to re-consider this "terminating" event as "true", the
            % frame needs to be shifted back 1. (if terminating event is not "true")
            % i.e.: [1, 1, 1, NaN; 47, 1, 2, 1] becomes
            %       [1, 1, 1, NaN; 46, 2, 1, NaN]
            if ~isnan(decompTrack(iDecompTrack).seqOfEvents(2,4))
                decompTrack(iDecompTrack).seqOfEvents(2,1) = decompTrack(iDecompTrack).seqOfEvents(2,1) - 1;
                decompTrack(iDecompTrack).seqOfEvents(1,2) = 1; % starting/initiating tracklet event
                decompTrack(iDecompTrack).seqOfEvents(2,2) = 2; % ending/terminating tracklet event
                decompTrack(iDecompTrack).seqOfEvents(2,4) = NaN; % true ending/terminating tracklet event
            end
            % name segment as "1"
            decompTrack(iDecompTrack).seqOfEvents(:,3) = ones(2,1);
            % Make sure all col4 of .seqOfEvents are NaN indicating "true"events
            if ~isnan(decompTrack(iDecompTrack).seqOfEvents(1,4))
                decompTrack(iDecompTrack).seqOfEvents(1,4) = NaN; % true ending/terminating tracklet event
            end
            
            %% To get tracksFeatIndxCG & tracksCoordAmpCG
            %%% Determine time frame of .seqOfEvents
            if totalInt == size(compTrack(iCompTrack).tracksFeatIndxCG,2)
                iniFrTracklet = decompTrack(iDecompTrack).seqOfEvents(1,1) - iniFrame + 1;
                endFrTracklet = decompTrack(iDecompTrack).seqOfEvents(2,1) - iniFrame + 1;
                % When a tracklet end by transition into another tracklet, we
                % no longer see detection
            else
                error("ERROR: tracksFeatIndxCG does not describe seqOfEvent correctly in compTrack!");
            end
            
            %%% Copy truncated tracksFeatIndxCG & tracksCoordAmpCG (to appropriate frame)
            % of each tracklet inside comptrack into decompTrack structure
            decompTrack(iDecompTrack).tracksFeatIndxCG = ...
                compTrack(iCompTrack).tracksFeatIndxCG(iTracklet,iniFrTracklet:endFrTracklet);
            decompTrack(iDecompTrack).tracksCoordAmpCG = ...
                compTrack(iCompTrack).tracksCoordAmpCG(iTracklet,((iniFrTracklet-1)*8+1):(endFrTracklet*8));
            % 20200925 TN: add handling to .aggregState
            if aggregField
                decompTrack(iDecompTrack).aggregState = ...
                    compTrack(iCompTrack).aggregState(iTracklet,iniFrTracklet:endFrTracklet);
            end
            
            %%% increment iDecompTrack after filling in each row of decompTrack structure
            iDecompTrack = iDecompTrack+1;
            
        end %(for iTracklet = 1:numTracklet)
        
    end %(if numTracklet == 1 ... else ...)
    
end %(for iCompTrack = 1:numCompTrack)

end