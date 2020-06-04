function EEG = filter_envelope(EEG, freq_lowcut, freq_highcut, order, method)
% filter_envelope(EEG, freq_lowcut, freq_highcut, order, method) - filter dataset.
%
% Usage:
%   >> [EEG] = filter_envelope(EEG, freq_cut);
%   >> [EEG] = filter_envelope(EEG, 10, 50, [], 'squared');
%   >> [EEG] = filter_envelope(EEG, [], [], [], 'hilbert');
%
% Inputs:
%   EEG          - input data struct
%
% Optional inputs:
%   freq_lowcut  - lowcut-off frequency (Hz). Not used for the hilbert
%                  transformation.
%   freq_highcut - order of the filter. Not used for the notch filter.
%                  1 element performs highpass, 2 elements perform bandpass
%   order        - order of the filter. Not used for the notch filter.
%                  The default is 3.
%   method       - type of filter ("abs","squared","hilbert")
%
% Outputs:
%   EEG          - output filtered data struct
%
% Note: this function uses filtfilt, butter, and hilbert from matlab
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

% Check inputs
if nargin < 2
    freq_lowcut = 10;
    freq_highcut = 50;
    order = 3;
    method = 'squared';
end

if nargin < 5
    order = 3;
    method = 'squared';
end

% Check envelope method
if ~any(strcmp(method,{'abs','squared','hilbert'}))
    method = 'squared';
end
disp(['Envelope method: ' method])

% Filter data
for ds = 1:numel(EEG)
    % Nyquist freq
    Fnyq  = EEG{ds}.srate/2;
    
    % Inver signal to column matrix
    if size(EEG{ds}.data,1) < size(EEG{ds}.data,2)
        EEG{ds}.data = EEG{ds}.data';
    end
    
    % Highpass/Bandpass data
    if ~isempty(freq_highcut)
        if numel(freq_highcut) == 1
            [b_filt,a_filt] = butter(order, freq_highcut/Fnyq,'high');
        elseif numel(freq_highcut) == 2
            [b_filt,a_filt] = butter(order, freq_highcut/Fnyq,'bandpass');
        else
            error('ERROR: freq_highcut can have either 1 or 2 elements!')
        end
        data_filt = detrend(filtfilt(b_filt, a_filt, EEG{ds}.data));
    else
        data_filt = detrend(EEG{ds}.data);
    end
    
    % Rectification
    if strcmp(method,'squared')
        % Squaring for rectifing: gain of 2 for maintianing the same energy 
        % in the output
        data_filt = 2*data_filt.*data_filt;
    elseif strcmp(method,'abs')
        data_filt = abs(data_filt);
    elseif strcmp(method,'hilbert')
        data_filt = abs(hilbert(data_filt));
    end
    
    % Low Pass
    if ~strcmp(method,'hilbert')
        [b_filt, a_filt] = butter(order,freq_lowcut/Fnyq,'low');
        data_filt = filtfilt(b_filt, a_filt, data_filt);
        
        if strcmp(method,'squared')
            data_filt = abs(sqrt(data_filt));
        end
    end
    
    % Band Pass
    EEG{ds}.data = data_filt';
end

% EOF