function [EEG] = get_epochs_tf_mtft(EEG, event_name, offset_before, offset_after, ...
                                    win_size, win_step, tapers, pad, freq_range, ...
                                    norm_factor, unit, chosen_trials)
% get_epochs_tf(EEG, event_name,offset_before,offset_after,chosen_sessions)
% - divides the EEG data by epochs and convert them in time-frequency using
%   multitaper from chronux_2_12.
%
% Usage:
%   >> EEG = get_epochs_tf_mtft(EEG, event_name, offset_before, offset_after, win_size, win_step);
%
% Inputs:
%   EEG            - input data struct. It contains EEG.data (channels x samples)
%   event_name     - name of the event that divide the epochs (str)
%   offset_before  - time before the event. It is in seconds (float)
%   offset_after   - time after the event. It is in seconds (float)
%   win_size       - Size of the sliding window for computing the stft.
%                    It is in seconds (float)
%   win_step       - Step of the sliding window for computing the stft.
%                    It is in seconds (float)
%
% Optional inputs:
%   tapers      - precalculated tapers from dpss or in the one of the
%                 following forms:
%                 (1) A numeric vector [TW K] where TW is the
%                     time-bandwidth product and K is the number of
%                     tapers to be used (less than or equal to
%                     2TW-1).
%                 (2) A numeric vector [W T p] where W is the
%                     bandwidth, T is the duration of the data and p
%                     is an integer such that 2TW-p tapers are used. In
%                     this form there is no default i.e. to specify
%                     the bandwidth, you have to specify T and p as
%                     well. Note that the units of W and T have to be
%                     consistent: if W is in Hz, T must be in seconds
%                     and vice versa. Note that these units must also
%                     be consistent with the units of params.Fs: W can
%                     be in Hz if and only if params.Fs is in Hz.
%                 The default is to use form 1 with TW=4 and K=7.
%                 Note that T has to be equal to win_size.
%	pad		    - (padding factor for the FFT) - optional (can take values
%                 -1,0,1,2...). -1 corresponds to no padding, 0 corresponds
%                 to padding to the next highest power of 2 etc.
%                 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0,
%                 to 512 points, if pad=1, we pad to 1024 points etc.
%			      we pad the FFT. The defaults to 0.
%   freq_range     - Min and max frequencies of the spectogram. ([float float])
%                    The default is all the range
%   norm_factor    - normalization factor for each channel. (n_channels)
%                    The default is [].
%   unit           - data unit as "power" or "db". The default is "power"
%   chosen_trials  - choosen sessions to concatinate data from. The default
%                    is all the sessions
%
% Outputs:
%   EEG            - output data struct with data divided by epochs.
%                    EEG.epochs (channels x samples x frequency x trials)
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

if nargin < 6
    error('ERROR: not enough inputs!');
end
if nargin < 7 || isempty(tapers)
    tapers = [4, 7];
end
if nargin < 8 || isempty(pad)
    pad = 0;
end
if nargin < 9 || isempty(freq_range)
    freq_range = [];
end
if nargin < 10 || isempty(norm_factor)
    norm_factor = [];
end
if nargin < 11 || isempty(unit)
    unit = 'power';
end
if nargin < 12
    chosen_trials = 1:numel(EEG);
end

for ds = chosen_trials
    % Invert data to row matrix
    EEG{ds}.data = transpose_data(EEG{ds}.data,'row');
    n_channels = size(EEG{ds}.data,1);
    n_samples = size(EEG{ds}.data,2);
    
    % Convert before and after epoch to samples
    offset_before_smpl = offset_before * EEG{ds}.srate;
    offset_after_smpl = offset_after * EEG{ds}.srate;
    % Convert windows to samples
    win_size_smpl = round(win_size * EEG{ds}.srate);
    win_step_smpl = round(win_step * EEG{ds}.srate);
    
    
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
        % Frequency info
        NFFT = max([256, 2^nextpow2(win_size_smpl)]);
        df = EEG{ds}.srate/NFFT;
        sfreqs = 0:df:EEG{ds}.srate/2; % all possible frequencies
        if ~isempty(freq_range)
            freq_idx = sfreqs>=freq_range(1) & sfreqs<=freq_range(2);
            f = sfreqs(freq_idx);
        else
            freq_range = [sfreqs(1), sfreqs(end)];
            f = sfreqs;
        end
        
        % Time info
        t = -offset_before_smpl:win_step_smpl:offset_after_smpl-1;
        t = t/EEG{ds}.srate;
        
        % epoch matrix size estimation and preallocation
        n_f = numel(f);         % calculate the number of unique fft points
        n_l = numel(t);         % calculate the number of signal frames
        epochs = nan(n_channels,n_l,n_f,n_events); % preallocate the epoch matrix
        % Select events samples
        events_smpl = cell2mat({EEG{ds}.event(events).latency});
        events_all_smpl = cell2mat({EEG{ds}.event.latency});
        
        event_counter = 0;
        for iEv = 1:numel(events_smpl)
            % Update data epoch
            data_interval = events_smpl(iEv) - offset_before_smpl - round(win_size_smpl/2) : events_smpl(iEv) + offset_after_smpl + round(win_size_smpl/2) - 1;
            [MTFT,~,~] = mtft(EEG{ds}.data(:,data_interval),EEG{ds}.srate,win_size_smpl,win_step_smpl,tapers,pad,freq_range,norm_factor,unit);
            epochs(:,:,:,iEv) = MTFT;
            
            find_events_in_epoch = ismember(events_all_smpl,events_smpl(iEv) - offset_before_smpl:events_smpl(iEv) + offset_after_smpl);
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
        EEG{ds}.times = t;
    else
        disp(['Trial ' num2str(ds) ' does NOT contain events within the data range...'])
    end
end

% EOF