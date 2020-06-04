function [STFT,f,t] = stft(data,srate,win_size, win_step, ...
                           freq_range, NFFT, window, ...
                           norm_factor, unit)
% short_time_ft(data,srate,win_size, win_step, freq_range) - apply short
% time Fourier transform
%
% Usage:
%   >> [STFT,f,t] = short_time_ft(data,srate,win_size, win_step);
%   >> [STFT,f,t] = short_time_ft(data,srate,win_size, win_step, freq_range, NFFT, window);
%
% Inputs:
%   data        - input data (n_channels x n_samples)
%   srate       - sampling frequency
%   win_size    - Size of the sliding window for computing the pmtm. In samples
%   win_step    - Step of the sliding window for computing the pmtm. In samples.
%
% Optional inputs:
%   freq_range  - Min and max frequencies of the spectogram. ([float float])
%                 The default is all the range.
%   NFFT        - Length of the signal for the FFT analisys. The default 
%                 is max(256, 2**nextpow2(n_samples)).
%   window      - spectral window. The default is the hamming window.
%   norm_factor - normalization factor for each channel. (n_channels).
%                 The default is [].
%   unit        - select unit as "power" or "db". The default is "power".
% 
% Outputs:
%   STFT        - spectrogram (n_channels x n_samples x n_freq)
%   f           - frequency vector (Hz)
%   t           - time vector (seconds)
%
% Note: this function uses stft, hamming, pow2db from matlab
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    error('ERROR: not enough inputs!'); 
end
if nargin < 5 || isempty(freq_range)
    freq_range = [];
end
if nargin < 6 || isempty(NFFT)
    NFFT = [];
    window = [];
end
if nargin < 7 || isempty(window)
    window = hamming(win_size,'periodic')./norm(hamming(win_size,'periodic'));
else
    window = window./norm(window);
end
if nargin < 8 || isempty(norm_factor)
    norm_factor = [];
end
if nargin < 9
    unit = 'power';
end

% Transpose data to column vector
data = transpose_data(data,'column');

n_channels = size(data,2);
n_samples = size(data,1);

if isempty(NFFT)
    NFFT = max([256, 2^nextpow2(win_size)]);
end

% Frequency info
df = srate/NFFT;
sfreqs = 0:df:srate/2; % all possible frequencies
if ~isempty(freq_range)
    freq_idx = sfreqs>=freq_range(1) & sfreqs<=freq_range(2);
    f = sfreqs(freq_idx);
else
    f = sfreqs;
    freq_idx = 1:numel(f);
end

% Time info
t = win_size/2:win_step:n_samples-1-win_size/2;
t = t/srate;

% stft matrix size estimation and preallocation
n_f = numel(f);         % calculate the number of unique fft points
n_l = numel(t);         % calculate the number of signal frames
STFT = zeros(n_channels, n_f, n_l); % preallocate the stft matrix

if ~isempty(norm_factor) && (length(norm_factor) ~= n_f)
    error('Normalization factor must have the same frequency dimension as the STFT!')
else
    norm_factor = transpose_data(norm_factor,'column');
end

% STFT (via time-localized FFT)
for l = 0:n_l-1
    % windowing
    data_win = bsxfun(@times, data(1+l*win_step : win_size+l*win_step,:), window);
    
    % FFT
    X = fft(data_win, NFFT);
    pow  = abs(X) ./ sqrt(srate);
    
    % update of the stft matrix
    STFT(:, :, 1+l) = pow(freq_idx,:)';
end

STFT = permute(STFT,[1,3,2]);

% Correct by normalization
if ~isempty(norm_factor)
    STFT = STFT./permute(repmat(norm_factor,[1,n_l,n_channels]),[3,2,1]);
end

if n_channels == 1
    STFT = squeeze(STFT);
end

% Convert amplitude spectrum to dB
if strcmpi(unit,'db')
    STFT = pow2db(STFT);
end

% EOF