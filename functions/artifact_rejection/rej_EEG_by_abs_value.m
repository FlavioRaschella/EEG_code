function [EEG] = rej_EEG_by_abs_value(EEG, rej_value, n_chan)
% rej_epochs_abs_value(EEG, rej_value) - set EEG data with absolute value
% higher than rej_value as "bad data".
%
% Usage:
%   >> EEG = rej_EEG_by_abs_value(EEG);
%   >> EEG = rej_EEG_by_abs_value(EEG, rej_value, n_chan);
%
% Inputs:
%   EEG       - input data struct. It contains EEG.data (channels x samples)
% 
% Optional inputs:
%   rej_value - threshold for selecting data as bad data.
%               The default is 200 (uV).
%   n_chan    - number of channels showing the same values over rej_value
%               in order to set these values as artifacts.
%               The default is 1.
%
% Outputs:
%   EEG       - output data struct with bad_data save in EEG.rej.rej_abs_val
%               (channels x samples)
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

if nargin < 2 || isempty(rej_value)
    rej_value = 200;
end
if nargin < 3 || isempty(n_chan)
    n_chan = 1;
end

for ds = 1:numel(EEG)
    % Invert data to row matrix
    data_tmp = transpose_data(EEG{ds}.data,'column');
    n_samples  = size(data_tmp,1);
    
    bad_data_count = sum(abs(data_tmp) > rej_value,2);
    bad_data       = bad_data_count >= n_chan;
    % Remove the signal around as well
    junk_len = EEG{ds}.srate; % 1s junk
    if sum(bad_data) ~= 0
        bad_data_start = find(diff(bad_data)==1)+1;
        bad_data_stop = find(diff(bad_data)==-1);
        
        for iBd = 1:numel(bad_data_start)
            if bad_data_start(iBd)-junk_len>0
                bad_data(bad_data_start(iBd)-junk_len:bad_data_start(iBd)) = 1;
            else
                bad_data(1:bad_data_start(iBd)) = 1;
            end
        end
        for iBd = 1:numel(bad_data_stop)
            if bad_data_stop(iBd)+junk_len<n_samples
                bad_data(bad_data_stop(iBd):bad_data_stop(iBd)+junk_len) = 1;
            else
                bad_data(bad_data_stop(iBd):end) = 1;
            end
        end
    else
        disp(['Trial ' num2str(ds) ' does NOT contain data over the set threshold...'])
    end
    
    % Update rejection structure
    if ~isfield(EEG{ds},'reject') || isempty(EEG{ds}.reject)
        EEG{ds}.reject = struct();
    end
    EEG{ds}.reject.rej_abs_val = bad_data;
end

% EOF