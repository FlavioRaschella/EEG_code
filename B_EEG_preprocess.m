% This script shows the pipeline for data processing.
% 0. Choose and load files
% 1. EOG artifact regression
% 2. Channel Interpolation
% 3. Filter data
% 4. Re-reference signal
% 5. Identify bad data
% 6. Separate data by trials
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean workspace
clear all; close all; clc

%% 0. Choose and load files
path_name = [pwd filesep 'data'];
EEG_files = {'EEG_data_prepared.mat'};
data_num = numel(EEG_files);

EEG = cell(data_num,1);
for fl = 1:data_num
    tmp = load([path_name filesep EEG_files{fl}]);
    EEG{fl} = tmp;
end

%% 1. EOG artifact regression
EEG = remove_EOG(EEG);

%% 2. Channel Interpolation

% To get Triangulation Weight :
% Cartool merge EEG + xzy > click xyz window > options> Exp triangulation distances
[~,~,~,tri_Mat] = readfile_sef([path_name filesep 'TriangulationWeights.sef']);
EEG = interpolate_channel(EEG,tri_Mat);

%% 3. Filter data
% bandpass (1-200Hz) and notch (50Hz)
EEG = filter_data(EEG,'band',[1, 200]);
EEG = filter_data(EEG,'notch',50);

%% Plot data
% eeglab must be on the path
% pop_eegplot(EEG{1}, 1, 1, 1);

%% 4. Re-reference signal
EEG = CAR(EEG);

%% 5. Identify bad data
EEG = rej_EEG_by_abs_value(EEG);

%% 6. Separate data by trials
EEG_rest  = get_trials(EEG, '"rest"');
EEG_trial = get_trials(EEG, '"trial"');

%% 7. Save data
save([path_name filesep 'EEG_data_processed.mat'],'EEG_rest','EEG_trial')

% EOF