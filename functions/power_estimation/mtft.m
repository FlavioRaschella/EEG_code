function [S,f,t] = mtft(data,srate,win_size,win_step,tapers,pad,freq_range,norm_factor,unit)
% mtft(data,srate,win_size,win_step,tapers,pad,freq_range,norm_factor,unit)
% - apply Multi-taper time-frequency spectrum
%
% Usage:
%   >> [S,f,t] = mtft(data,srate,win_size,win_step,params);
%   >> [S,f,t] = mtft(data,srate,win_size,win_step,params,norm_factor,unit);
%
% Inputs:
%   data        - input data (n_channels x n_samples)
%   srate       - sampling frequency
%   win_size    - Size of the sliding window for computing the pmtm. In samples
%   win_step    - Step of the sliding window for computing the pmtm. In samples.
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
%   freq_range  - Min and max frequencies of the spectogram. ([float float])
%                 The default is 0 and srate/2.
%   NFFT        - Length of the signal for the FFT analisys. The default
%                 is max(256, 2**nextpow2(n_samples)).
%   window      - spectral window. The default is the hamming window.
%   norm_factor - normalization factor for each channel. (n_channels).
%                 The default is [].
%   unit        - select unit as "power" or "db". The default is "power".
%
% Outputs:
%   S           - spectrogram (n_channels x n_samples x n_freq)
%   f           - frequency vector (Hz)
%   t           - time vector (seconds)
%
% Note: this function is adapted from chronux_2_12
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    error('ERROR: not enough inputs!');
end
if nargin < 5 || isempty(tapers)
    tapers = [4, 7];
end
if nargin < 6 || isempty(pad)
    pad = 0;
end
if nargin < 7 || isempty(freq_range)
    freq_range = [];
end
if nargin < 8 || isempty(norm_factor)
    norm_factor = [];
end
if nargin < 9
    unit = 'power';
end

% Set params
params = struct();
params.Fs = srate;
if ~isempty(tapers)
    params.tapers = tapers;
end
if ~isempty(freq_range)
    params.fpass = freq_range;
end
if ~isempty(pad)
    params.pad = pad;
end

[tapers,pad,Fs,fpass,~,~,params]=getparams(params);
if length(params.tapers)==3 && win_size/Fs~=params.tapers(2)
    error('ERROR: duration of data in params.tapers is inconsistent with win_size, modify params.tapers(2) to proceed!')
end

% Transpose data to column vector
data     = transpose_data(data,'column');
[n_samples,n_channels]   = size(data);
win_size = round(win_size); % number of samples in window
win_step = round(win_step); % number of samples to step through
nfft     = max([256, 2^(nextpow2(win_size)+pad)]);
f        = getfgrid(Fs,nfft,fpass);
params.tapers = dpsschk(tapers,win_size,Fs); % check tapers

% Set time windows
t = win_size/2:win_step:n_samples-1-win_size/2;
t = t/srate;

% mtft matrix preallocation
n_l = length(t);
n_f = length(f);
S = zeros(n_channels,n_f,n_l);

if ~isempty(norm_factor) && (length(norm_factor) ~= n_f)
    error('Normalization factor must have the same frequency dimension as the STFT!')
else
    norm_factor = transpose_data(norm_factor,'column');
end

% MTFT
for l = 0:n_l-1
    data_win = data(1+l*win_step : win_size+l*win_step,:);
    [s,f] = mtspectrumc(data_win,params);
    S(:,:,1+l) = s';
end

S = permute(S,[1,3,2]);

% Correct by normalization
if ~isempty(norm_factor)
    S = S./permute(repmat(norm_factor,[1,n_l,n_channels]),[3,2,1]);
end

if n_channels == 1
    S = squeeze(S);
end

% Convert amplitude spectrum to dB
if strcmpi(unit,'db')
    S = pow2db(S);
end

% EOF