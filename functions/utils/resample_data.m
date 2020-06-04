function EEG = resample_data(EEG,Fs_resample,fc,df)
% resample_data() - resample dataset.
%
% Usage:
%   >> [EEG] = resample_data( EEG, Fs_resample);
%   >> [EEG] = resample_data( EEG, Fs_resample, fc, df);
%
% Inputs:
%   EEG         - input data struct. EEG.data can also contain epochs.
%                 (channels x samples) or (channels x samples x epochs)
%   Fs_resample - frequency to resample (Hz)
%
% Optional inputs:
%   fc          - anti-aliasing filter cutoff (pi rad / sample). The default
%                 is 0.9.
%   df          - anti-aliasing filter transition band width (pi rad /
%                 sample). The default is 0.2
%
% Outputs:
%   EEG        - output data struct with resampled EEG.data.
%
% Note: this function uses resample from matlab and pop_firwsord, firws
% from eeglab
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('ERROR: data and its frequency should be inputed.')
end

% Default cutoff frequency (pi rad / smp)
if nargin < 3 || isempty(fc)
    fc = 0.9;
end

% Default transition band width (pi rad / smp)
if nargin < 4 || isempty(df)
    df = 0.2;
end

% Check that fc is in range
if fc < 0 || fc > 1
    error('ERROR: Anti-aliasing filter cutoff freqeuncy out of range (<0 || >1).')
end

if ~iscell(EEG)
    EEG = {EEG};
end

for ds = 1:numel(EEG)
    % finding the best ratio
    [p,q] = rat(Fs_resample/EEG{ds}.srate, 1e-12);
    n_channels = size(EEG{ds}.data,1);
    n_samples = size(EEG{ds}.data,2);
    FS_old = EEG{ds}.srate;
        
    fprintf('resampling dataset %d: %3.4f Hz\n', ds, EEG{ds}.srate*p/q);
    
    EEG{ds}.srate   = EEG{ds}.srate*p/q;
    if numel(size(EEG{ds}.data)) == 2
        EEG{ds}.data = resample_padding( double(transpose_data(EEG{ds}.data,'column')), p, q, fc, df )';
    elseif numel(size(EEG{ds}.data)) == 3
        data_tmp = nan(n_channels,EEG{ds}.srate*n_samples/FS_old,size(EEG{ds}.data,3));
        for iEp = 1:size(EEG{ds}.data,3)
            data_tmp(:,:,iEp) = resample_padding( double(transpose_data(EEG{ds}.data(:,:,iEp),'column')), p, q, fc, df )';
        end
        EEG{ds}.data = data_tmp;
    else
        error('ERROR: EEG.data has dimensions ~= 2 or 3!')
    end
    EEG{ds}.pnts  = size(EEG{ds}.data,2);
    EEG{ds}.xmax  = EEG{ds}.xmin + (EEG{ds}.pnts-1)/EEG{ds}.srate;
    EEG{ds}.times = linspace(EEG{ds}.xmin, EEG{ds}.xmax, EEG{ds}.pnts);
    try EEG{ds}.epoch_before = EEG{ds}.epoch_before*EEG{ds}.srate/FS_old; ,catch ,end
    try EEG{ds}.epoch_after  = EEG{ds}.epoch_after*EEG{ds}.srate/FS_old; ,catch ,end
    
    % resample events
    if numel(size(EEG{ds}.data)) == 2
        if isfield(EEG{ds}, 'event') && isfield(EEG{ds}.event, 'type') && ischar(EEG{ds}.event(1).type)
            latency_resampled = round(cell2mat({EEG{ds}.event.latency})*(EEG{ds}.srate/FS_old));
            for iEl = 1:numel(latency_resampled)
                EEG{ds}.event(iEl).latency = latency_resampled(iEl);
            end
        else
            disp('no events to resample.')
        end
    else
        if isfield(EEG{ds}, 'epoch') && isfield(EEG{ds}.epoch, 'eventlatency') && ~isempty(EEG{ds}.epoch(1).eventlatency)
            for iEp = 1:numel(EEG{ds}.epoch)
                for iEl = 1:numel(EEG{ds}.epoch(iEp).eventlatency)
                    EEG{ds}.epoch(iEp).eventlatency{iEl} = round(EEG{ds}.epoch(iEp).eventlatency{iEl}*(EEG{ds}.srate/FS_old));
                end
            end
        end
    end
end

% store dataset
disp('resampling finished/n');

return;

function data_resample = resample_padding(data, p, q, fc, df)

if length(data) < 2
    data_resample = data;
    return;
end

% padding to avoid artifacts at the beginning and at the end
% The resample_data command introduces substantial artifacts at beginning
% and end of data when raw data show DC offset (e.g. as in DC recorded
% continuous files) when MATLAB Signal Processing Toolbox is present (and
% MATLAB resample.m command is used).
% Even if this artifact is short, it is a filtered DC offset and will be
% carried into data, e.g. by later highpass filtering to a substantial
% amount (easily up to several seconds).
% The problem can be solved by padding the data at beginning and end by a
% DC constant before resampling.

% Conservative custom anti-aliasing FIR filter
nyq = 1 / max([p q]);
fc = fc * nyq; % Anti-aliasing filter cutoff frequency
df = df * nyq; % Anti-aliasing filter transition band width
m = pop_firwsord('kaiser', 2, df, 0.002); % Anti-aliasing filter kernel
b = firws(m, fc, windows('kaiser', m + 1, 5)); % Anti-aliasing filter kernel
b = p * b; % Normalize filter kernel to inserted zeros

% Padding
nPad = ceil((m / 2) / q) * q; % Datapoints to pad, round to integer multiple of q for unpadding
startPad = repmat(data(1, :), [nPad 1]);
endPad = repmat(data(end, :), [nPad 1]);

% Resampling
data_resample = resample([startPad; data; endPad], p, q, b);

% Remove padding
nPad = nPad * p / q; % # datapoints to unpad
data_resample = data_resample(nPad + 1:end - nPad, :); % Remove padded data