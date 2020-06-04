function [EEG_trial] = get_trials(EEG, event_start, event_stop)
% get_trials(EEG, event_start) - divide the EEG data by trials.
%
% Usage:
%   >> [EEG_trial] = get_trials(EEG, event_start);
%
% Inputs:
%   EEG         - input data struct
%   event_start - name of the event that divide the trials (start trial)
% 
% Optional inputs:
%   event_stop  - name of the events that stop the trials (stop trial).
%                 The default is empty.
%
% Outputs:
%   EEG_trial   - output data struct with data divided by trials
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

if nargin < 3
    event_stop = [];
end

EEG_trial = [];

for ds = 1:numel(EEG)
    % Find events
    events_st = strcmpi({EEG{ds}.event.type},event_start);
    if ~isempty(event_stop)
        events_end = strcmpi({EEG{ds}.event.type},event_stop);
        if (length(events_st) ~= length(events_end))
            error('ERROR: event_stop and event_start have different length!');
        end
        if any(events_end-events_st<0)
            error('ERROR: some event_stop are smaller than event_start!');
        end
    end
    
    % Divide trials
    if ~isempty(events_st)
        events_st_smpl = cell2mat({EEG{ds}.event(events_st).latency});
        if ~isempty(event_stop)
            events_end_smpl = cell2mat({EEG{ds}.event(events_end).latency});
        else
            events_end_smpl = [events_st_smpl(2:end), EEG{ds}.pnts];
        end
        
        for iTr = 1:numel(events_st_smpl)
            EEG_tmp = EEG{ds};
            % Update data
            EEG_tmp.data = EEG_tmp.data(:,events_st_smpl(iTr):events_end_smpl(iTr));
            % Update trial
            EEG_tmp.trials = [EEG_tmp.trials, iTr];
            % Update points, time, xmax
            EEG_tmp.pnts = size(EEG_tmp.data,2);
            EEG_tmp.times = 1/EEG_tmp.srate:1/EEG_tmp.srate:length(EEG_tmp.data)/EEG_tmp.srate;
            EEG_tmp.xmax = length(EEG_tmp.data)/EEG_tmp.srate;
            
            % Update EOG data
            EEG_tmp.EOG = EEG_tmp.EOG(:,events_st_smpl(iTr):events_end_smpl(iTr));
            
            % Reset event structure
            EEG_tmp.event = struct('type',[],...
                'latency',[],...
                'urevent',[]);
            
            if isfield(EEG_tmp,'reject') && isfield(EEG_tmp.reject,'rej_abs_val') ...
                    && length(EEG_tmp.reject.rej_abs_val) == length(EEG{ds}.data)
                EEG_tmp.reject.rej_abs_val = EEG_tmp.reject.rej_abs_val(events_st_smpl(iTr):events_end_smpl(iTr));
            end
            
            % Store EEG_tmp
            EEG_trial = [EEG_trial, {EEG_tmp}];
        end
    end
end

% EOF