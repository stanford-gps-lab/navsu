function jd = epochs2jd(epochs)
% epochs2jd  Converts epoch to Julian date. Non-vectorized version.
% Version: 28 Sep 03
% Usage:   jd=gps2jd(epochs)
% Input:   epochs  - GPS epoch
% Output:  jd      - Julian date

% Copyright (c) 2011, Michael R. Craymer
% All rights reserved.
% Email: mike@craymer.com

jdgps = navsu.time.cal2jd(1980,1,6);             % beginning of GPS week numbering
jd = jdgps + epochs/3600/24;