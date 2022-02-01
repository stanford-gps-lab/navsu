function jd = cal2jd_vect(yr, mn, dy)
% CAL2JD  Converts calendar date to Julian date using algorithm
%   from "Practical Ephemeris Calculations" by Oliver Montenbruck
%   (Springer-Verlag, 1989). Uses astronomical year for B.C. dates
%   (2 BC = -1 yr). Non-vectorized version. See also DOY2JD, GPS2JD,
%   JD2CAL, JD2DOW, JD2DOY, JD2GPS, JD2YR, YR2JD.
% Version: 2011-11-13
% Usage:   jd=cal2jd(yr,mn,dy)
% Input:   yr - calendar year (4-digit including century)
%          mn - calendar month
%          dy - calendar day (including factional day)
% Output:  jd - jJulian date

% Copyright (c) 2011, Michael R. Craymer
% All rights reserved.
% Email: mike@craymer.com

if nargin ~= 3
  warning('Incorrect number of input arguments');
  return;
end

% make all inputs are row vectors of equal sizes
nTimes = max([length(yr) length(mn) length(dy)]);
colVec = any([size(yr, 1)>1, size(mn, 1)>1, size(dy,  1)>1]);

dy = makeRowVec(dy, nTimes);
mn = makeRowVec(mn, nTimes);
yr = makeRowVec(yr, nTimes);

invalidEpochs = mn < 1 ...
              | mn > 12 ...
              | (mn == 2 & dy > 29) ...
              | (ismember(mn, [4 6 9 11]) & dy > 30) ...
              | dy > 31;

date1 = 4.5 + 31*(10 + 12*1582);   % Last day of Julian calendar (1582.10.04 Noon)
date2 = 15.5 + 31*(10 + 12*1582);  % First day of Gregorian calendar (1582.10.15 Noon)
date = dy + 31*(mn + 12*yr);

i = mn <= 2;
yr(i) = yr(i) - 1;
mn(i) = mn(i) + 12;

% number of leap years?
b = NaN(size(date));
b(date <= date1) = -2;
i = date >= date2; 
b(i) = fix(yr(i)/400) - fix(yr(i)/100);
% if any(isnan(b))
%   warning('cal2jd_vect: Dates between October 5 & 15, 1582 do not exist');
% %   return;
% end

jd = fix(365.25*yr) + fix(30.6001*(mn+1)) + b + 1720996.5 + dy;
i = yr <= 0;
jd(i) = fix(365.25*yr(i)-0.75) + fix(30.6001*(mn(i)+1)) + b(i) + 1720996.5 + dy(i);

% account for invalid inputs
jd(invalidEpochs) = NaN;

if colVec
    % transform back
    jd = jd';
end
end

function x = makeRowVec(x, n)
% Turns input into row vector of length n.
if numel(x) == 1
    x = repmat(x, 1, n);
elseif iscolumn(x)
    x = x';
end

end
