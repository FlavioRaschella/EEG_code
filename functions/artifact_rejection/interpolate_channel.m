function [EEG] = interpolate_channel(EEG,tri_Mat,threshold,f_noise)
% interpolate_channel(EEG,tri_Mat,threshold,f_noise) - find and interpolate
% noisy channels
%
% Usage:
%   >> [EEG, ix_noise ,noisy_ch] = interpolate_channel(EEG,tri_Mat,threshold,f_noise);
%   >> [EEG, ix_noise ,noisy_ch] = interpolate_channel(EEG,tri_Mat);
%
% Inputs:
%   EEG            - input data struct
%   tri_Mat        - Triangulation matrix to identify electrodes positions. 
%
% Optional inputs:
%   threshold      - Multiplier of the signal's poswer std to define it a signal
%                    as noisy. The default is 3.
%   f_noise        - frequency above which identify noisy signals. The default 
%                    is 30Hz.
%
% Outputs:
%   EEG            - output intepolated data struct
%
% Note: this function uses nextpow2, pwelch from Matlab
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(EEG)
    EEG = {EEG};
end

if nargin < 4
    threshold = 3;
end

for ds = 1:numel(EEG)
    % Invert data to column matrix
    data_tmp = transpose_data(EEG{ds}.data,'column');
    srate = EEG{ds}.srate;
    % Number of channels
    n_channels = size(data_tmp,2);
    
    % pwelch window
    win_welch = 2^(nextpow2(srate));
    % frequency axis
    f_axis = 0:srate/win_welch:srate/2;
    % power spectral density estimate (PSD)
    PSD = pwelch(data_tmp,win_welch);
    
    f_max = floor(f_axis(rank(cov(PSD.'))));
    
    % Look at power above 30Hz to identify noisy channnels
    if nargin < 5
        % in the case of low passed signals, use f_max/2 if this value is 
        % lower than 30Hz.
        f_noise = min([f_max/2,30]);
    end
    
    PSD_noise = mean(PSD(f_axis>=f_noise,:));
    [PSD_noise, PSD_chan] = sort(PSD_noise,'descend');
    PSD_noise_diff = diff(PSD_noise);
    
    % PSD noise distribution
    PSD_noise_std = PSD_noise-median(PSD_noise);
    SD = quantile(abs(PSD_noise_std),0.6827);
    PSD_noise_std = PSD_noise_std./SD;
    
    % Get noisy channels
    PSD_chan_noise = PSD_chan(PSD_noise_std > threshold);
    
    % Plot channels' power
    figure();
    subplot(2,1,1); 
    plot(f_axis,PSD);
    hold on; plot(f_axis,PSD(:,PSD_chan_noise),'r.')
    [~, ix_alpha] = min(abs(f_axis-10));
    axis([1, 120, 0, 2*max(PSD(ix_alpha,:))]);
    ylabel('PSD magnitude')
    xlabel('Frequency')
    title('Welch Power Sprectal Density Estimate')
    subplot(2,1,2); 
    plot(1:size(PSD,2),PSD_noise);
    hold on; plot(1:size(PSD_chan_noise,2),PSD_noise(1:size(PSD_chan_noise,2)),'r')
    hold on; plot(1:size(PSD_chan_noise,2),PSD_noise(1:size(PSD_chan_noise,2)),'rx')
    ylabel('Magnitude of noise in power')
    xlabel('Channels (descending order by noise)')
    set(gcf,'Position',[0,20,400,980])
    
    % Display identified noisy channels
    noisy_ch = PSD_chan_noise;
    disp([num2str(numel(noisy_ch)) ' noisy channels were detected; ' num2str(round(100*numel(noisy_ch)/n_channels,1)) '% of the dataset.'])
    
    % Interpolate channels that were detected as noisy
    tri_Mat(diag(true(size(tri_Mat,1),1))) = 0;
    interp_mat = tri_Mat;
    valid_ch = true(size(data_tmp,2),1);
    valid_ch(PSD_chan_noise) = false;
    interp_mat(~valid_ch,:) = false;
    
    bad_ch = false(size(data_tmp,2),1);
    n_neighbors = sum(logical( interp_mat(:,PSD_chan_noise)));
    
    % If there is no valid neighbor for a channel, this is obmitted in the 
    % first interpolation step and interpolated in a second step after 
    % interpolating the neighbors
    
    bad_ch(PSD_chan_noise(n_neighbors==0)) = true;
    PSD_chan_noise(n_neighbors==0) = [];
    n_neighbors(n_neighbors==0) = [];
    
    if any(PSD_chan_noise)
        % Interpolation
        data_tmp(:,PSD_chan_noise) = data_tmp *  interp_mat(:,PSD_chan_noise) ./ n_neighbors;
        
        % Try to fix bad channels
        interp_mat = tri_Mat;
        interp_mat(bad_ch,:) = false;
        n_neighbors = sum(logical(interp_mat(:,bad_ch)));
        if any(bad_ch)
            if any(n_neighbors)
                bad_ix = find(bad_ch);
                data_tmp(:,bad_ch) = data_tmp * interp_mat(:,bad_ix(n_neighbors>0)) ./ n_neighbors(n_neighbors>0);
                bad_ch(bad_ix(n_neighbors>0)) = false;
                if any(bad_ch)
                    fprintf([num2str(sum(bad_ch)) ' channels could not be interpolated due to no valid neighbor electrodes \n' ...
                        'setting these channels as "bad" '])
                    data_tmp(:,bad_ch) = 0;
                end
            else
                fprintf([num2str(sum(bad_ch)) ' channels could not be interpolated due to no valid neighbor electrodes \n' ...
                    'setting these channels as "bad" '])
                data_tmp(:,bad_ch) = 0;
            end
        end
    end
    
    % Store data
    EEG{ds}.data = data_tmp';
    EEG{ds}.noisy_ch = noisy_ch;
    EEG{ds}.bad_ch = bad_ch(bad_ch~=0);
end

% EOF