function [zmap_nc,zmap_cluster,zmap_min_max] = permutation_test(tf,n_permutes,p_val,plot_data,frex,times2save)
% permutation_test(EEG, data_ch) -sStatistics via permutation testing on
% time-frequency data.
%
% Usage:
%   >> [zmap_nc,zmap_cluster,zmap_min_max] = permutation_test(tf);
%   >> [zmap_nc,zmap_cluster,zmap_min_max] = permutation_test(tf,n_permutes,pval,plot_data,frex,times2save);
%
% Inputs:
%   tf         - input time-frequency array (samples x frequencies x trials)
%
% Optional inputs:
%   n_permutes - number of permutations. The default is 1000.
%   p_val      - p-value for t/z-test. The default is 0.05.
%   plot_data  - set whether to plot the result images. The default is false.
%   frex       - frequencies for plotting data. The default is 1 : numel(frequencies).
%   times2save - time for plotting data. The default is 1 : numel(samples).
%
% Outputs:
%   EEG        - output filtered data struct
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(numel(size(tf)) == 4)
assert(size(tf,1) == 2)

num_frex = size(tf,2);
num_time = size(tf,3);
ntrials = size(tf,4);

if nargin<2 || isempty(n_permutes)
    n_permutes = 1000;
end

if nargin<3 || isempty(p_val)
    p_val = 0.05;
end

if nargin<4 || isempty(plot_data)
    plot_data = False;
end
if nargin<5 || isempty(frex)
    frex = 1:num_frex;
end
if nargin<6 || isempty(times2save)
    times2save = 1:num_time;
end

clim = [0 20000];
xlim = [-.1 1]; % for plotting

% convert p-value to Z value
zval = abs(norminv(p_val));

% initialize null hypothesis maps
permmaps = zeros(n_permutes,num_frex,num_time);

% for convenience, tf power maps are concatenated
%   in this matrix, trials 1:ntrials are from channel "1"
%   and trials ntrials+1:end are from channel "2"
diffmap = squeeze(mean(tf(2,:,:,:),4 )) - squeeze(mean(tf(1,:,:,:),4 ));
tf3d = cat(3,squeeze(tf(1,:,:,:)),squeeze(tf(2,:,:,:)));


% generate maps under the null hypothesis
for permi = 1:n_permutes
    
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(size(tf3d,3));
    temp_tf3d = tf3d(:,:,randorder);
    
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = squeeze( mean(temp_tf3d(:,:,1:ntrials),3) - mean(temp_tf3d(:,:,ntrials+1:end),3) );
end

%% Show thresholded maps non-corrected for multiple comparison

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;
zmap_nc = zmap;

% now some plotting...
if plot_data
    figure
    
    subplot(221)
    imagesc(times2save,frex,diffmap);
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','nor')
    title('TF map of real power values')
    
    subplot(222)
    imagesc(times2save,frex,diffmap);
    hold on
    contour(times2save,frex,logical(zmap),1,'linecolor','k');
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
    title('Power values and outlined significance regions')
    
    subplot(223)
    imagesc(times2save,frex,zmap);
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    set(gca,'clim',[-10 10],'xlim',xlim,'ydir','no')
    title('Thresholded TF map of Z-values')
    
    suptitle('non-corrected map')
    
end
%% Corrections for multiple comparisons

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
    
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
end

%% Show histograph of maximum cluster sizes

% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*p_val));

if plot_data
    figure
    hist(max_cluster_sizes,20);
    hold on; plot([cluster_thresh,cluster_thresh],[0,max(hist(max_cluster_sizes,20))]);
    xlabel('Maximum cluster sizes'), ylabel('Number of observations')
    title('Expected cluster sizes under the null hypothesis')
end

%% Plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end

zmap_cluster = zmap;

if plot_data
    % plot tresholded results
    figure
    subplot(221)
    imagesc(times2save,frex,diffmap)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('TF power, no thresholding')
    set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
    
    subplot(222)
    imagesc(times2save,frex,diffmap)
    hold on
    contour(times2save,frex,logical(zmap_cluster),1,'linecolor','k')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('TF power with contour')
    set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','norm')
    
    subplot(223)
    imagesc(times2save,frex,zmap_cluster)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('z-map, thresholded')
    set(gca,'clim',[-13 13],'xlim',xlim,'ydir','normal')
    
    suptitle('cluster-corrected map')
end
%% Now with max-pixel-based thresholding

% find the threshold for lower and upper values
thresh_lo = prctile(max_val(:,1),    100*p_val); % what is the
thresh_hi = prctile(max_val(:,2),100-100*p_val); % true p-value?
% note about the above code: a 2-tailed test actually requires pval/2 on each tail;
% thus, the above is actually testing for p<.1 !

% threshold real data
zmap = diffmap;
zmap(zmap>thresh_lo & zmap<thresh_hi) = 0;

zmap_min_max = zmap;

if plot_data
    figure
    subplot(221)
    imagesc(times2save,frex,diffmap)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('tf power map, no thresholding')
    set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','n')
    
    subplot(222)
    imagesc(times2save,frex,diffmap)
    hold on
    contour(times2save,frex,logical(zmap_min_max),1,'linecolor','k')
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('tf power map with contour')
    set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','normal')
    
    subplot(223)
    imagesc(times2save,frex,zmap_min_max)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('tf power map, thresholded')
    set(gca,'clim',[-mean(clim) mean(clim)],'xlim',xlim,'ydir','no')
    
    suptitle('min-max corrected map')
end

% EOF