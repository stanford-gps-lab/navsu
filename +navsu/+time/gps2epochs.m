function epochs = gps2epochs(gpsWeek,tow)
% gps2epochs
% DESCRIPTION:
%   Convert from  GPS week and GPS time of week to GPS epochs (seconds 
%   since start of GPS time)
% INPUT:
%   gpsWeek  = GPS week with no rollover [Nx1]
%   gpsTow   = GPS time of week [Nx1]
%    -- or --
%   [gpsWeek gpsTow] = same as above, only in a single matrix of size [Nx2]
%
% OUTPUT:
%   epochs = GPS epoch [Nx1]
%
% See also: navsu.time.epochs2gps

if nargin == 1 && size(gpsWeek, 2) == 2, % input is [weekNum TOW] (size N x 2)
    tow = gpsWeek(:, 2);
    gpsWeek = gpsWeek(:, 1);
end

epochs = gpsWeek*604800+tow;

end