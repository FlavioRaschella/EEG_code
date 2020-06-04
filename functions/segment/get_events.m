function [EEG] = get_events(EEG, events, events_name, events_multiplier)
% get_trials(EEG, events, events_name) - get events as store them in the EEG struct.
%
% Usage:
%   >> EEG = get_events(EEG, events);
%
% Inputs:
%   EEG               - input data struct. No cell.
%   events            - indices of the events (n_event)
%   events_name       - name of the events that divide the trials
% 
% Optional inputs:
%   events_multiplier - events multiplier. It is used to convert the events
%                       in the sample timeframe of the data. The default is
%                       empty --> not used. (float)
%
% Outputs:
%   EEG               - output data struct with assigned events
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(EEG)
    error('ERROR: EEG must be a struct! You inputed a cell.')
end

if nargin < 4
    events_multiplier = [];
end

if ~isfield(EEG,'event') || ~isfield(EEG.event,'urevent') ...
        || ~isfield(EEG.event,'type') || ~isfield(EEG.event,'latency')
    EEG.event(length(events)) = struct('type',[],...
        'latency',[],...
        'urevent',[]);
    event_start_idx = 0;
elseif (isfield(EEG,'event') && isfield(EEG.event,'urevent') ...
        && isfield(EEG.event,'type') && isfield(EEG.event,'latency')) ...
        && (isempty(EEG.event(1).urevent) || isempty(EEG.event(1).type) || ...
        isempty(EEG.event(1).latency))
    event_start_idx = 0;
else
    event_start_idx = length(EEG.event);
end
    
% Store events
for i = 1 : length(events)
    EEG.event(event_start_idx+i).urevent = event_start_idx+i;
    EEG.event(event_start_idx+i).type    = events_name;
    if isempty(events_multiplier)
        EEG.event(event_start_idx+i).latency = round(events(i));
    else
        EEG.event(event_start_idx+i).latency = round(events(i)*events_multiplier);
    end
end

% Sort events
EEG = sort_events(EEG);

% EOF