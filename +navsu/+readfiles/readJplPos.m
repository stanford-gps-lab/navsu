function Peph = readJplPos(filename, startRow, endRow)
% readJplClock
% DESCRIPTION:
%   Parser for JPL real-time orbit products. 
%
% INPUT:
%   filename    - the name of the file (please include the path as well)
%
% OUTPUT:
%   Peph        - structure containing parsed data- built to match
%                 navsu.readfiles.readSp3- check that for more info.
%
% See also: navsu.readfiles.readJplClock, navsu.readfiles.readSp3

%% Initialize variables.
delimiter = {' ','S'};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: categorical (%C)
%	column2: categorical (%C)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%C%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
data = table(dataArray{1:end-1}, 'VariableNames', {'E','constellation',...
    'svn','epoch_jpl','event','pos_x','pos_y','pos_z','vel_x','vel_y','vel_z','E05','E4','E5','E09','E6'});


%% convert and rearrange various things
% convert from jpl epochs to gps epochs
epochs = navsu.time.jpl2epochs(data.epoch_jpl);

% constellations
% this doesn't really do anything right now
constellation = ones(size(epochs));

% convert from svn to prn
prns = navsu.svprn.svn2prn(data.svn,epochs,constellation(1));

% some data is not available at all here
clock_bias  = nan(size(epochs));
clock_drift = nan(size(epochs));

% I think this is the event flag?
Event = data.event;

position = [data.pos_x data.pos_y data.pos_z]*1000;
velocity = [data.vel_x data.vel_y data.vel_z]*1000;


%% Put into Peph struct format
Peph.PRN           = prns;
Peph.clock_bias    = clock_bias;
Peph.clock_drift   = clock_drift;
Peph.position      = position;
Peph.velocity      = velocity;
Peph.Event         = Event;
Peph.epochs        = epochs;
Peph.constellation = constellation;















