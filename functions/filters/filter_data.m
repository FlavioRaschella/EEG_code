function EEG = filter_data(EEG, type, freq_cut, order)
% filter_BP(EEG, type, freq_cut, order) - filter dataset.
%
% Usage:
%   >> [EEG] = filter_data( EEG, type, freq_cut);
%   >> [EEG] = filter_data( EEG, type, freq_cut, order);
%
% Inputs:
%   EEG      - input data struct
%   type     - type of filter ("band","high","low","notch")
%   freq_cut - cut-off frequency (Hz)
%
% Optional inputs:
%   order    - order of the filter. Not used for the notch filter. 
%              The default is 3.
%
% Outputs:
%   EEG      - output filtered data struct
%
% Note: this function uses filtfilt, butter, and iirnotch from matlab
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

% Check inputs
if nargin < 4
    order = 3;
end

if ~any(strcmp(type,{'band','high','low','notch'}))
    error('ERROR: type can only be "band","high","low", or "notch"!')
end

for ds = 1:numel(EEG)
    % Nyquist freq
    Fnyq  = EEG{ds}.srate/2;
    
    switch type
        case 'band'
            if numel(freq_cut) ~= 2
                error('ERROR: freq_cut for bandpass filter must have 2 elements! (e.g. [1 40])')
            end
            [b_filt,a_filt] = butter(order, freq_cut/Fnyq,'bandpass');
        case 'high'
            if numel(freq_cut) ~= 1
                error('ERROR: freq_cut for highpass filter must have 1 elements! (e.g. [1])')
            end
            [b_filt,a_filt] = butter(order, freq_cut/Fnyq,'high');
        case 'low'
            if numel(freq_cut) ~= 1
                error('ERROR: freq_cut for lowpass filter must have 1 elements! (e.g. [100])')
            end
            [b_filt,a_filt] = butter(order, freq_cut/Fnyq,'low');
        case 'notch'
            if numel(freq_cut) ~= 1
                error('ERROR: freq_cut for notch filter must have 1 elements! (e.g. [50])')
            end
            % Design a filter with a Q-factor of Q=20 to remove a freq_cut Hz tone
            [b_filt,a_filt] = iirnotch(freq_cut/Fnyq, freq_cut/Fnyq/20);
            
            % Aternative computation of Notch filter
            % d  = fdesign.notch('N,F0,BW',4,freq_cut/Fnyq,1/Fnyq);
            % Hd = design(d);
            % [b_filt,a_filt] = sos2tf(Hd.sosMatrix);
            
        otherwise
            error('ERROR: switch type different from the possible choises!')
    end
    
    % Invert data to column matrix
    EEG{ds}.data = transpose_data(EEG{ds}.data,'column');
    
    % Band Pass
    EEG{ds}.data = filtfilt(b_filt, a_filt, EEG{ds}.data)';    
end

% EOF