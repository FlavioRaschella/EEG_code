function EEG = remove_EOG(EEG, EOG_ch)
% remove_EOG(EEG, data_ch) - remove the EOG artefacts using regression.
%
% Usage:
%   >> [EEG] = remove_EOG(EEG);
%   >> [EEG] = remove_EOG(EEG, data_ch);
%
% Inputs:
%   EEG      - input data struct. EEG.data (channels x signals)
%
% Optional inputs:
%   EOG_ch   - channels to use as EOG if the EOG recording are missing (channels). 
%              The default is [].
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

% Select channels to use for re-referencing
if nargin<2
    EOG_ch = [];
end

for ds = 1:numel(EEG)    
    % Invert data to column matrix
    EEG{ds}.data = transpose_data(EEG{ds}.data,'column');
    
    if isfield(EEG{ds},'EOG') && ~isempty(EEG{ds}.EOG)
        EOG_tmp = transpose_data(EEG{ds}.EOG,'column');
        % Check data and EOG length
        if size(EEG{ds}.data,1) ~= size(EOG_tmp,1)
            error('ERROR: EOG and data have different length!')
        end
    elseif ~isempty(EOG_ch)
        EOG_tmp = EEG{ds}.data(:,EOG_ch);
    else
        error('ERROR: either EOG should be in the EEG structure or some EEG channels sould be used as EOG!')
    end
    
    % EEG regression with EOG
    EEG{ds}.data = EEG{ds}.data - EOG_tmp*(EOG_tmp\EEG{ds}.data);
end

% EOF