function data = transpose_data(data,direction)
% transpose_data(data,direction) - transpose the data in the selected direction.
%
% Usage:
%   >> data = transpose_data(data);
%   >> data = transpose_data(data,direction);
%
% Inputs:
%   data      - input data matrix.
%
% Optional inputs:
%   direction - transpose the data in 'row' or 'column'. The default is
%               'column'.
%
% Outputs:
%   data      - output transposed data
%
% author of this script: Flavio Raschella
% last modify: 25/05/2020
% contact: flavio.raschella@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    direction = 'column';
end
if ~any(strcmpi(direction,{'row','column'}))
    error('ERROR: direction can either be "row" or "column"!')
end

% Transpose
if strcmpi(direction,'column') && (size(data,1) < size(data,2))
    data = data';
end

if strcmpi(direction,'row') && (size(data,1) > size(data,2))
    data = data';
end

% EOF