function EEG = zero_mean(EEG, data_ch)
% zero_mean(EEG,data_ch) - re-reference the data their own mean.
%
% Usage:
%   >> [EEG] = zero_mean(EEG);
%   >> [EEG] = CAR(EEG, data_ch);
%
% Inputs:
%   EEG      - input data struct. EEG.data (channels x signals)
% 
% Optional inputs:
%   data_ch  - channels to re-reference (channels). The default is all data
%              channels
% 
% Outputs:
%   EEG      - output filtered data struct
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

for ds = 1:numel(EEG)    
    % Invert data to column matrix
    EEG{ds}.data = transpose_data(EEG{ds}.data,'column');
    
    % Select channels to use for re-referencing
    if nargin<2 || isempty(data_ch)
        data_ch = 1:size(EEG{ds}.data,2);
    end
    
    % Zero-mean
    ZM = mean(EEG{ds}.data(:,data_ch),1);
    
    % Re-referenced data to CA
    EEG{ds}.data = bsxfun(@minus, EEG{ds}.data(:,data_ch), ZM)';
    
    % Save info
    EEG{ds}.ref = 'zero_mean';
    EEG{ds}.ref_signal = ZM;
    EEG{ds}.ref_chan = data_ch;
end

% EOF