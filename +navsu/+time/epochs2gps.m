function [gpsWeek, gpsTow] = epochs2gps(epochs)
% epochs2gps 
% DESCRIPTION:
%   Convert from GPS epochs (seconds since start of GPS time) to GPS week
%   and GPS time of week
% INPUT:
%   epochs = GPS epoch [Nx1]
%
% OUTPUT:
%   gpsWeek  = GPS week with no rollover [Nx1]
%   gpsTow   = GPS time of week [Nx1]
%
% See also: navsu.time.gps2epochs

gpsWeek = floor(epochs./(7*86400));
gpsTow = epochs-gpsWeek*7*86400;


end