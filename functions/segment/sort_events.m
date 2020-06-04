function [EEG] = sort_events(EEG)
% order_events(EEG) - sort events in ascending order.
%
% Usage:
%   >> EEG = order_events(EEG);
%
% Inputs:
%   EEG               - input data struct.
%
% Outputs:
%   EEG               - output data struct with sorted events
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
    
% Store events
for ds = 1 : length(EEG)
    [~, sort_idx] = sort(cell2mat({EEG{ds}.event.latency}));
    latency_tmp   = cell2mat({EEG{ds}.event(sort_idx).latency});
    type_tmp      = {EEG{ds}.event(sort_idx).type};
    % Check latencies ascending order
    if any(diff(latency_tmp)<0)
        error('ERROR: issue with sorting the latencies. Not in ascending order!')
    end
    % Store events
    for iEl = 1:numel(sort_idx)
        EEG{ds}.event(iEl).latency = latency_tmp(iEl);
        EEG{ds}.event(iEl).type    = type_tmp{iEl};
        EEG{ds}.event(iEl).urevent = iEl;
    end
end

% Adjust if not cell
if ~iscell_flag
    EEG = EEG{1};
end

% EOF