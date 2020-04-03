function jd = epochs2jd(epochs)

% GPS2JD  Converts epoch to Julian date. Non-vectorized version.
%   See also CAL2JD, DOY2JD, JD2CAL, JD2DOW, JD2DOY, JD2GPS,
%   JD2YR, YR2JD.
% Version: 28 Sep 03
% Usage:   jd=gps2jd(gpsweek,sow,rollover)
% Input:   gpsweek  - GPS week number
%          sow      - seconds of week since 0 hr, Sun (default=0)
%          rollover - number of GPS week rollovers (default=0)
% Output:  jd       - Julian date

% Copyright (c) 2011, Michael R. Craymer
% All rights reserved.
% Email: mike@craymer.com

jdgps = navsu.time.cal2jd(1980,1,6);             % beginning of GPS week numbering
jd = jdgps + epochs/3600/24;