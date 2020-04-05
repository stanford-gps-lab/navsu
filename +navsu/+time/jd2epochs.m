function epochs = jd2epochs(jds)
% jd2epochs
% DESCRIPTION: 
%   Convert from julian date to GPS epochs (seconds since start of GPS
%   time)
% INPUT:
%   jds = Julian date vector [Nx1]
%
% OUTPUT:
%   epochs = GPS epochs vector [Nx1]
%
% See also: navsu.time.epochs2jd

[weeks,tow] = navsu.time.jd2gps(jds);

epochs = weeks*604800+tow;



end