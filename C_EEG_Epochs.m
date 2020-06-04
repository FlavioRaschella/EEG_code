% This script shows the pipeline for separating the data in epochs, both.
% 0. Reading preprocessed EEG files
% 1. Load events for each trial
% 2. Extract epochs for one event
% 3. Resample data
% 4. Extract epochs for one event as time-frequency
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean workspace
clear all; close all; clc

%% 0. Reading preprocessed EEG files
path_name = [pwd filesep 'data'];
EEG_files = 'EEG_data_processed.mat';
load([path_name filesep EEG_files])

%% 1. Load events for each trial
event_file = {'trial1_events.mat';...
              'trial2_events.mat'};
          
assert(numel(event_file) == numel(EEG_trial))

for iFl = 1:numel(event_file)
    event_tmp = load([path_name filesep event_file{iFl}]);
    
    EEG_trial{iFl} = get_events(EEG_trial{iFl}, event_tmp.event1.start, 'event1_start');
    EEG_trial{iFl} = get_events(EEG_trial{iFl}, event_tmp.event1.stop,  'event1_stop');
    EEG_trial{iFl} = get_events(EEG_trial{iFl}, event_tmp.event2.start, 'event2_start');
    EEG_trial{iFl} = get_events(EEG_trial{iFl}, event_tmp.event2.stop,  'event2_stop');
end

%% 2. Extract epochs for one event
EEG_epo = get_epochs(EEG_trial, 'event1_start', 1, 2);

% Choose whether normalise the data
% EEG_epo_norm = epochs_normalize(EEG_epo);

% Extract epoch stats
EEG_epo = get_epochs_stats(EEG_epo);

%% 3. Resample data
EEG_epo_res = resample_data(EEG_epo, 128);

%% 4. Extract epochs for one event as time-frequency
EEG_epo_tf = get_epochs_tf_stft(EEG_trial, 'event1_start', 1, 2, 0.5, 0.002);
EEG_epo_tf = get_epochs_stats(EEG_epo_tf);

% EOF