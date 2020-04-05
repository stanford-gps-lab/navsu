function mjd = jd2mjd(jd)
% jd2mjd
% DESCRIPTION:
%   Convert from julian date to modified julian date
% INPUT:
%   jd = Julian date [Nx1]
%
% OUTPUT:
%   mjd = modified Julian date [Nx1]
%
% See also: navsu.time.mjd2jd

% convert from modified julian date to julian date
mjd= jd-2400000.5;

end