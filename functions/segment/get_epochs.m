function [EEG] = get_epochs(EEG, event_name, offset_before, offset_after, chosen_trials)
% get_epochs(EEG, event_name,offset_before,offset_after,chosen_sessions) - 
% divides the EEG data by epochs.
%
% Usage:
%   >> EEG = get_epochs(EEG, event_name, offset_before, offset_after);
%   >> EEG = get_epochs(EEG, event_name, offset_before, offset_after, chosen_sessions);
%
% Inputs:
%   EEG            - input data struct. It contains EEG.data (channels x samples)
%   event_name     - name of the event that divide the epochs (str)
%   offset_before  - time before the event. It is in seconds (int)
%   offset_after   - time after the event. It is in seconds (int)
% 
% Optional inputs:
%   chosen_trials  - choosen sessions to concatinate data from. The default
%                    is all the sessions.
%
% Outputs:
%   EEG            - output data struct with data divided by epochs.
%                    EEG.epochs (channels x samples x trials)
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

if nargin < 5
    chosen_trials = 1:numel(EEG);
end

for ds = chosen_trials
    % Invert data to row matrix
    EEG{ds}.data = transpose_data(EEG{ds}.data,'row');
    n_channels = size(EEG{ds}.data,1);
    n_samples = length(EEG{ds}.data);
    
    % Convert before and after epoch to samples
    offset_before_smpl = round(offset_before * EEG{ds}.srate);
    offset_after_smpl  = round(offset_after * EEG{ds}.srate);
    
    % Update epoch structure
    EEG{ds}.epoch = struct('event',[]);
    event_fields = fieldnames(EEG{ds}.event);
    for iFld = 1:numel(event_fields)
        EEG{ds}.epoch.(['event' event_fields{iFld}]) = [];
    end
    
    % Find events
    events_tmp = find(strcmpi({EEG{ds}.event.type},event_name));
    % Check that epochs around the events are within the data range, and
    % that data in epochs are good
    if ~isempty(events_tmp)
        events = [];
        for iEv = 1:numel(events_tmp)
            if (EEG{ds}.event(events_tmp(iEv)).latency - offset_before_smpl > 0) && ...
                    (EEG{ds}.event(events_tmp(iEv)).latency + offset_after_smpl < n_samples)
                if isfield(EEG{ds},'reject') && isfield(EEG{ds}.reject,'rej_abs_val') ...
                        && length(EEG{ds}.reject.rej_abs_val) == n_samples
                    epoch_interval = EEG{ds}.event(events_tmp(iEv)).latency - offset_before_smpl:...
                        EEG{ds}.event(events_tmp(iEv)).latency + offset_after_smpl-1;
                    if ~any(EEG{ds}.reject.rej_abs_val(epoch_interval))
                        events = [events, events_tmp(iEv)];
                    else
                        disp(['Trial ' num2str(ds) '; event ' events_tmp(iEv) ' is in a bad epoch...'])
                    end
                else
                    events = [events, events_tmp(iEv)];
                end
            end
        end
        n_events = length(events);
    else
        disp(['Trial ' num2str(ds) ' does NOT contain the event: '  event_name '...'])
        continue
    end
    
    % Divide trials
    if ~isempty(events)
        % Assign epoch array
        epoch_len = offset_before_smpl+offset_after_smpl;
        epochs = nan(n_channels,epoch_len,n_events);
        % Select events samples
        events_smpl = cell2mat({EEG{ds}.event(events).latency});
        events_all_smpl = cell2mat({EEG{ds}.event.latency});
        
        event_counter = 0;
        for iEv = 1:numel(events_smpl)
            % Update data epoch
            epoch_interval = events_smpl(iEv) - offset_before_smpl + 1 : events_smpl(iEv) + offset_after_smpl;
            epochs(:,:,iEv) = EEG{ds}.data(:,epoch_interval);
            
            find_events_in_epoch = ismember(events_all_smpl,epoch_interval);
            
            % Update epoch struct
            EEG{ds}.epoch(iEv).event = event_counter + 1: event_counter + sum(find_events_in_epoch);
            
            for iFld = 1:numel(event_fields)
                if strcmpi(event_fields{iFld},'latency')
                    EEG{ds}.epoch(iEv).(['event' event_fields{iFld}]) = num2cell((cell2mat({EEG{ds}.event(find_events_in_epoch).(event_fields{iFld})}) - events_smpl(iEv))/EEG{ds}.srate);
                else
                    EEG{ds}.epoch(iEv).(['event' event_fields{iFld}]) = {EEG{ds}.event(find_events_in_epoch).(event_fields{iFld})};
                end
            end
            
            % Update counter
            event_counter = event_counter + sum(find_events_in_epoch);
        end
        % Store epochs
        EEG{ds}.data = epochs;
        EEG{ds}.epoch_before = offset_before_smpl;
        EEG{ds}.epoch_after  = offset_after_smpl;
        EEG{ds}.xmin = -offset_before_smpl/EEG{ds}.srate;
        EEG{ds}.xmax = (offset_after_smpl-1)/EEG{ds}.srate;
        EEG{ds}.times = EEG{ds}.xmin:1/EEG{ds}.srate:EEG{ds}.xmax;
    else
        disp(['Trial ' num2str(ds) ' does NOT contain events within the data range...'])
    end
end

% EOF