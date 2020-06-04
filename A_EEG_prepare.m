% This script loads the data and store them into an EEG structure.
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean workspace
clear all; close all; clc

%% Set data path info
% Folder Path
path_name = [pwd filesep 'data'];
% EEG data files
EEG_files = {'EEG_data.mat'};
% Get data num
data_num = numel(EEG_files);

%% Load and process data
EEG = cell(data_num,1);
for fl = 1:data_num
    eeg_tmp = load([path_name filesep EEG_files{fl}]);
    
    options = {'events',eeg_tmp.Trig_on,'events_name',eeg_tmp.Event,'eog',[1,3],'trial',fl};
    EEG{fl} = load_EEG_struct(path_name,EEG_files{fl},eeg_tmp.data,eeg_tmp.fs,options);
end
disp ('Data loaded!')
%% Plot data, check for abnormal signals

for fl = 1:data_num
    waterplot(EEG{fl},[],true)
end

%% Save data
for fl = 1:data_num
    eeg_data = EEG{fl};
    save([path_name filesep EEG_files{fl}(1:end-4) '_prepared.mat'],'-struct','eeg_data')
end

% EOF