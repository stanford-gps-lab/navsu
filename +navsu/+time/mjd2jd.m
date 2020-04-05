function jd = mjd2jd(mjd)
% mjd2jd
% DESCRIPTION:
%   Convert from modified julian date to julian date
% INPUT:
%   mjd = modified Julian date [Nx1]
%
% OUTPUT:
%   jd = Julian date [Nx1]
%
% See also: navsu.time.jd2mjd

% convert from modified julian date to julian date
jd = mjd+2400000.5;

end