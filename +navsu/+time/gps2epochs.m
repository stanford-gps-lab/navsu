function epochs = gps2epochs(gpsWeek,tow)
% epochs2gps 
% DESCRIPTION:
%   Convert from  GPS week and GPS time of week to GPS epochs (seconds 
%   since start of GPS time)
% INPUT:
%   gpsWeek  = GPS week with no rollover [Nx1]
%   gpsTow   = GPS time of week [Nx1]
%
% OUTPUT:
%   epochs = GPS epoch [Nx1]
%
% See also: navsu.time.epochs2gps

epochs = gpsWeek*604800+tow;

end