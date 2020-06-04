function [EEG] = epochs_normalize(EEG,offset_before,norm_factor,chosen_trials)
% epochs_normalize(EEG,offset_before) - normalize EEG epochs.
%
% Usage:
%   >> EEG = epochs_normalize(EEG);
%   >> EEG = epochs_normalize(EEG, offset_before, norm_factor, chosen_trials);
%
% Inputs:
%   EEG            - input data struct. It contains EEG.data that can either be
%                    (channels x samples x epochs) or
%                    (channels x samples x frequencies x epochs)
%
% Optional inputs:
%   offset_before  - time before the event to use for normalization. It is
%                    in seconds (int). The default is EEG.offset_before
%   norm_factor    - different normalization factor. 
%                    In case of 3D data, norm_factor should be (channels).
%                    In case of 3D data, norm_factor could either be 
%                    (frequencies) or (channels x frequencies).
%                    The default is [].
%   chosen_trials  - choosen sessions to concatenate data from. The default
%                    is all the sessions.
%
% Outputs:
%   EEG            - output data struct with normalized epochs.
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

if nargin < 2
    offset_before = [];
end
if nargin < 3
    norm_factor = cell(1,numel(EEG));
end
if nargin < 4
    chosen_trials = 1:numel(EEG);
end

if ~isempty(norm_factor) && numel(norm_factor)~=numel(EEG)
    error('ERROR: norm_factor has different length from EEG!')
end

for ds = chosen_trials
    data_size  = size(EEG{ds}.data);
    data_dim   = numel(data_size);
    n_channels = data_size(1);
    n_samples  = data_size(2);
    n_epochs   = data_size(end);
    
    % Get normalization factor
    if isempty(norm_factor{ds})
        % Get offset_before for baseline
        if isempty(offset_before)
            offset_before_smpl = EEG{ds}.epoch_before;
        else
            offset_before_smpl = offset_before * EEG{ds}.srate;
        end
        if data_dim == 3
            norm_factor{ds} = repmat(mean(EEG{ds}.data(:,1:offset_before_smpl,:),2),[1,n_samples,1]);
        elseif data_dim == 4
            norm_factor{ds} = repmat(mean(EEG{ds}.data(:,1:offset_before_smpl,:,:),2),[1,n_samples,1,1]);
        end
    else
        if data_dim == 3 && numel(norm_factor{ds}) ~= n_channels
            error('ERROR: norm_factor must have the same lenght as n_channels!');
        elseif data_dim == 4 && all(numel(norm_factor{ds}) ~= [data_size(3), data_size(1)*data_size(3)])
            error('ERROR: norm_factor must have the same lenght as n_frequency or size (n_channels x n_frequency)!');
        end
        if data_dim == 4 && numel(norm_factor{ds}) == data_size(3)
            norm_factor{ds} = permute(repmat(transpose_data(norm_factor{ds}),[1,n_samples,n_channels,n_epochs]),[3,2,1,4]);
        elseif data_dim == 4 && numel(norm_factor{ds}) == data_size(1)*data_size(3)
            norm_factor{ds} = permute(repmat(norm_factor{ds},[1,1,n_samples,n_epochs]),[1,3,2,4]);
        end
    end
    EEG{ds}.data = (EEG{ds}.data-norm_factor{ds})./norm_factor{ds};
end

% EOF