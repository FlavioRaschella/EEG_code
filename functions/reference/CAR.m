function EEG = CAR(EEG, data_ch, CAR_ch)
% CAR(EEG, data_ch, CAR_ch) - re-reference the data to a common average (CA).
%
% Usage:
%   >> [EEG] = CAR(EEG);
%   >> [EEG] = CAR(EEG, data_ch, CAR_ch);
%
% Inputs:
%   EEG      - input data struct. EEG.data (channels x signals)
%
% Optional inputs:
%   data_ch  - channels to re-reference (channels). The default is all data
%              channels
%   CAR_ch   - channel to use for re-referencing (channel). The default is 
%              all data channels
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
    if nargin<2
        data_ch = 1:size(EEG{ds}.data,2);
    end
    if nargin<3
        CAR_ch = data_ch;
    end
    if isempty(data_ch)
        data_ch = 1:size(EEG{ds}.data,2);
    end
    
    % Get CA
    CA = mean(EEG{ds}.data(:,CAR_ch),2);
    
    % Re-referenced data to CA
    EEG{ds}.data = bsxfun(@minus, EEG{ds}.data(:,data_ch), CA)';
    
    % Save info
    EEG{ds}.ref = 'CA';
    EEG{ds}.ref_signal = CA;
    EEG{ds}.ref_chan = data_ch;
end

% EOF