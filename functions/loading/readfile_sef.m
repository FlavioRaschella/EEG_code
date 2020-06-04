function [num_channels,num_time_frames,srate,data] = readfile_sef(filename)
% readfile_sef(filename) - reads the information contained in a .sef file
% from Cartool
%
% Usage:
%   >> [num_channels,num_time_frames,srate,data] = readfile_sef(filename);
%
% Inputs:
%   filename        - input data
%
% Outputs:
%   num_channels    - number of channels
%   num_time_frames - number of time frames
%   srate           - sampling rate
%   data            - array (num_time_frames x num_channels)
%
% Note: Cartool: http://brainmapping.unige.ch/Cartool.htm
% 
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open sef file
fid = fopen(filename,'r');

% Read fixed part of the header
version         = strcat(fread(fid,4,'int8=>char')');
num_channels    = fread(fid,1,'int32');
num_at_channels = fread(fid,1,'int32');
num_time_frames = fread(fid,1,'int32');
srate           = fread(fid,1,'float32');
t_year          = fread(fid,1,'int16');
t_month         = fread(fid,1,'int16');
t_day           = fread(fid,1,'int16');
t_hour          = fread(fid,1,'int16');
t_minute        = fread(fid,1,'int16');
t_second        = fread(fid,1,'int16');
t_millisecond   = fread(fid,1,'int16');

% Read channels
channels = strcat(fread(fid,[8,num_channels],'int8=>char')');

% Read data
data = fread(fid,[num_channels,num_time_frames],'float32')';

% Close file
fclose(fid);

% EOF