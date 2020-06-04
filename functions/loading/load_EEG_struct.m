function [EEG] = load_EEG_struct(path_name,file_name,data,srate,varargin)
% load_EEG_struct(path_name,file_name,data,srate,varargin) - load data into
% a struct.
%
% Usage:
%   >> [EEG] = load_EEG_struct(path_name,file_name,data,srate);
%   >> [EEG] = load_EEG_struct(path_name,file_name,data,srate,varargin);
%
% Inputs:
%   path_name   - name of the datafile path (str)
%   file_name   - name of the datafile (str)
%   data        - EEG data (n_channel x n_samples)
%   srate       - sampling rate (float)
%   varargin    - Optional inputs to assign to the EEG struct. 
%                 Format is : 'function', value
% 
% Optional inputs:
%   setname     - descriptive title for the dataset (str)
%   trial       - number of the trial (int)
%   eog         - EOG signals, or number of the EEG channels used as EOG
%                 ((n_channel x n_samples) or [int]). e.g. [1,2,3,4] 4 channels
%   events      - latency of the events ([float])
%   events_name - name of the events ([str])
%   events_type - set whether events are in 'time' or 'samples' (str)
%   ref         - reference (str)
%   ch_names    - name of the channels ({char} x n_channel)
%   ch_x        - channels x position ([float] x n_channel)
%   ch_y        - channels y position ([float] x n_channel)
%   ch_z        - channels z position ([float] x n_channel)
%
% Outputs:
%   EEG         - output EEG data struct
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that varargin contains an even number of elements
if iscell(varargin) && length(varargin) == 1
    varargin = varargin{1};
end
if rem(length(varargin),2)
    error('ERROR: varargin has an odd number of elemnets!')
end

setname = [];
trial = [];
eog = [];
events = [];
events_name = [];
events_type = 'samples';
ref = [];
ch_names = [];
ch_x = [];
ch_y = [];
ch_z = [];

% Loop along varargin
for i = 1:2:length(varargin)
    fn = varargin{i};
    val = varargin{i+1};
    
    switch lower(fn) % some fields have a special case
        case 'setname'
            setname = val;
        case 'trial'
            trial = val;
        case 'eog'
            eog = val;
        case 'events'
            events = val;
        case 'events_name'
            events_name = val;
        case 'ref'
            ref = val;
        case 'ch_names'
            ch_names = val;
        case 'ch_x'
            ch_x = val;
        case 'ch_y'
            ch_y = val;
        case 'ch_z'
            ch_z = val;
        otherwise
            disp(['Input ' fn ' not recognized! Skipped...'])
    end
end

%%%%%%%%%%%%%%%%%
% Store data info
%%%%%%%%%%%%%%%%%
EEG = struct();

% Path info
EEG.filepath = path_name;
EEG.filename = file_name;
disp(['loading EEG data from ' path_name filesep file_name '...'])

% Data
EEG.data = transpose_data(data,'row');
EEG.srate = srate;
EEG.nbchan = size(EEG.data,1);
EEG.pnts = size(EEG.data,2);
EEG.times = 1/EEG.srate:1/EEG.srate:length(EEG.data)/EEG.srate;

EEG.xmin = 0;
EEG.xmax = length(EEG.data)/EEG.srate;

%%%%%%%%%%%%%%%%%
% Optional values
%%%%%%%%%%%%%%%%%

% Setname
if ~isempty(setname)
    EEG.setname = setname;
else
    EEG.setname = 'Continuous data';
end

% Trial
if ~isempty(trial)
    EEG.trials = trial;
else
    EEG.trials = 0;
end

% Events
EEG.event = struct('type',[],...
    'latency',[],...
    'urevent',[]);
if ~isempty(events) && ~isempty(events_name) && (length(events) == length(events_name))
    for i = 1 : length(events)
        EEG.event(i).urevent = i;
        if iscell(events_name(i))
            EEG.event(i).type    = events_name{i};
        else
            EEG.event(i).type    = events_name(i);
        end
        if strcmpi(events_type,'samples')
            EEG.event(i).latency = events(i);
        elseif strcmpi(events_type,'time')
            EEG.event(i).latency = events(i)*EEG.srate;
        else
            error('ERROR: event_type can either be "samples" or "time"!')
        end
    end
    disp('events stored...')
else
    disp('no events stored...')
end

% EOG
if ~isempty(eog)
    if any(numel(eog) == [2,4])
        EEG.EOG = data(eog,:);
    elseif any(size(eog,2) == size(data))
        EEG.EOG = transpose_data(eog,'row');
    end
    disp('EOG stored...')
else
    EEG.EOG = [];
    disp('no EOG stored...')
end

% Reference
if ~isempty(ref)
    EEG.ref = ref;
else
    EEG.ref = 'common';
end

% Channels info
EEG.chanlocs(EEG.nbchan) = struct('theta',[],...
    'radius',[],...
    'labels',[],...
    'sph_theta',[],...
    'sph_phi',[],...
    'X',[],...
    'Y',[],...
    'Z',[],...
    'sph_radius',[],...
    'type',[],...
    'ref',[],...
    'urchan',[]);

for i = 1 : EEG.nbchan
    EEG.chanlocs(i).urchan = i;
end

if isempty(ch_names) || (numel(ch_names) == EEG.nbchan)
    ch_names = cellfun(@num2str,num2cell(1:EEG.nbchan),'un',0);
end
for i = 1 : EEG.nbchan
    try 
        EEG.chanlocs(i).labels = strtrim(ch_names{i});
    catch
        EEG.chanlocs(i).labels = strtrim(ch_names(i,:));
    end
end

if ~isempty(ch_x) && (numel(ch_x) == EEG.nbchan) && ...
   ~isempty(ch_y) && (numel(ch_y) == EEG.nbchan) && ...
   ~isempty(ch_z) && (numel(ch_z) == EEG.nbchan)
    for i = 1 : EEG.nbchan
        EEG.chanlocs(i).X = ch_x(i);
        EEG.chanlocs(i).Y = ch_y(i);
        EEG.chanlocs(i).Z = ch_z(i);
    end
end

disp(' ')

% EOF