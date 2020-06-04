function [EEG] = remove_events(EEG, events_name)
% remove_events(EEG, events_name) - remove events from the EEG struct.
%
% Usage:
%   >> EEG = remove_events(EEG, events_name);
%
% Inputs:
%   EEG               - input data struct. No cell.
%   events_name       - name of the events to remove (str)
% 
% Outputs:
%   EEG               - output data struct without removed events
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iscell_flag = 1;
if ~iscell(EEG)
    iscell_flag = 0;
    EEG = {EEG};
end
    
% Remove events
for ds = 1 : length(EEG)
    event_idx = strcmpi({EEG{ds}.event.type},events_name);
    
    if any(events_name)
        latency_tmp   = cell2mat({EEG{ds}.event(~event_idx).latency});
        type_tmp      = {EEG{ds}.event(~event_idx).type};
        % Store events
        for iEl = 1:numel(latency_tmp)
            EEG{ds}.event(iEl).latency = latency_tmp(iEl);
            EEG{ds}.event(iEl).type    = type_tmp{iEl};
            EEG{ds}.event(iEl).urevent = iEl;
        end
    end
end

% Adjust if not cell
if ~iscell_flag
    EEG = EEG{1};
end

% EOF