function [EEG] = get_epochs_stats(EEG, base_mean, base_std)
% get_epochs_stats(EEG, base_mean, base_std) - get stats on EEG epochs.
%
% Usage:
%   >> EEG = get_epochs_stats(EEG);
%   >> EEG = get_epochs_stats(EEG, base_mean, base_std);
%
% Inputs:
%   EEG       - input data struct. It contains EEG.data that can either be 
%               (channels x samples x epochs) or 
%               (channels x samples x frequencies x epochs)
%
% Optional inputs:
%   base_mean - baseline mean for normalization (channels) or 
%               (channels x frequencies). The default is [].
%   base_std  - baseline std for normalization (channels) or 
%               (channels x frequencies). The default is [].
% 
% Outputs:
%   EEG       - output data struct with epochs stasts:
%               EEG.epochstats.mean; EEG.epochstats.std; EEG.epochstats.sem; 
%               EEG.epochstats.n_epochs; EEG.epochstats.snr; 
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

if nargin < 2
    base_mean = cell(1,numel(EEG));
end
if nargin < 3
    base_std = cell(1,numel(EEG));
end

if numel(base_mean)~=numel(base_std)
    error('ERROR: normalization factors have different length!')
end
if numel(base_std)~=numel(EEG)
    error('ERROR: normalization factors have different length from EEG!')
end

for ds = 1:numel(EEG)
    % Get number of trials
    data_size  = size(EEG{ds}.data);
    data_dim   = numel(data_size);
    n_channels = data_size(1);
    n_samples  = data_size(2);
    n_epochs   = data_size(end);
    
    % Check data dimension
    if ~any(data_dim == [3,4])
        error('ERROR: data input dimension must be either 3 or 4!')
    end
    
    % Compute mean
    epoch_mean = mean(EEG{ds}.data,data_dim);
    % Compute std
    epoch_std = std(EEG{ds}.data,0,data_dim);
    % Compute sem
    epoch_sem = epoch_std/sqrt(n_epochs);
    
    % Create struct
    EEG{ds}.epochstats = struct('mean',epoch_mean,...
        'std',epoch_std,...
        'sem',epoch_sem,...
        'n_epochs',n_epochs);
    
    % Compute SNR
    if ~isempty(base_mean{ds}) && ~isempty(base_std{ds})
        if data_dim == 3
            if (length(base_mean{ds}) ~= n_channels) || (length(base_std{ds}) ~= n_channels)
                error('ERROR: base_mean and base_std must have the same length of the n_channels!')
            end
            epoch_snr = abs(epoch_mean - repmat(transpose_data(base_mean{ds}),[1,n_samples]))./...
                (epoch_std+repmat(transpose_data(base_std{ds}),[1,n_samples]));
        elseif data_dim == 4
            n_freq = data_size(3);
            if any(size(base_mean{ds}) ~= [n_channels,n_freq]) || any(size(base_std{ds}) ~= [n_channels,n_freq])
                error('ERROR: base_mean and base_std must have size (n_channels x n_freq)!')
            end
            epoch_snr = abs(epoch_mean - permute(repmat(transpose_data(base_mean{ds},'row'),[1,1,n_samples]),[1,3,2]))./...
                (epoch_std + permute(repmat(transpose_data(base_std{ds},'row'),[1,1,n_samples]),[1,3,2]));
        end
        EEG{ds}.epochstats.snr = epoch_snr;
    end
end

% EOF