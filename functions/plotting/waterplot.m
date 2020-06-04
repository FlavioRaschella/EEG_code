function waterplot(EEG, color, zscore_flag)
% waterplot(EEG, color, zscore_flag) - plot EEG data.
%
% Usage:
%   >> [EEG] = filter_data( EEG, type, freq_cut);
%   >> [EEG] = filter_data( EEG, type, freq_cut, order);
%
% Inputs:
%   EEG         - input data struct.
%
% Optional inputs:
%   color       - set a color for all the signals. The default is rainbow.
%   zscore_flag - set whether zscore the data. The default is false.
%
% Note: this function uses zscore from Matlab.
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(EEG)
    EEG = EEG{1};
end

% Transpose data to column vector
data = transpose_data(EEG.data,'column');
% Get number of channels
n_channels = size(data,2);

if nargin < 2
    color = repmat({'k';'b';'g';'m';'c';'r'},ceil(n_channels/6),1);
end

if nargin < 3
    zscore_flag = false;
end

if isempty(color)
    color = repmat({'k';'b';'g';'m';'c';'r'},ceil(n_channels/6),1);
elseif ischar(color)
    color = repmat({color},n_channels,1);
end

if zscore_flag == true
    data = zscore(data);
    detrend_flag = false;
else
    data = detrend(data);
    detrend_flag = true;
end

time_field = false;
if isfield(EEG,'times')
    time_field = true;
    t = EEG.times;
end

% Get data max
amplitude_max = round(median(max(data,[],1)));

figure(); hold all;
for ch = 1:n_channels
    if time_field
        plot(t,data(:,ch)+(ch-1)*amplitude_max*sign(amplitude_max),'Color',color{ch})
    else
        plot(data(:,ch)+(ch-1)*amplitude_max*sign(amplitude_max),'Color',color{ch})
    end
end

ylim([-amplitude_max*sign(amplitude_max) n_channels*amplitude_max*sign(amplitude_max)])

yticks(linspace(0,amplitude_max*(n_channels-1)*sign(amplitude_max),n_channels))
if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs(1).labels)
    yticklabels({EEG.chanlocs.labels})
else
    yticklabels(1:n_channels)
end

title(['Data; scale: ' num2str(amplitude_max) '; zscore: ' num2str(zscore_flag) '; detrend: ' num2str(detrend_flag)])
ylabel('Channels')
if time_field
    xlabel('Time')
else
    xlabel('Samples')
end

% EOF