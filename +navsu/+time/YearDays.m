function nDays = YearDays(year)
% yearDays
% DESCRIPTION:
%   Compute the number of days in a given year accounting for leap years
% INPUT:
%   year = the year [1x1]
%
% OUTPUT:
%   nDays = number of days in the year [1x1]
%
% See also: monthly calendars, daily calendars with jokes or comics

% In the Gregorian calendar 3 criteria must be taken into account to identify leap years:
% The year is evenly divisible by 4;
% If the year can be evenly divided by 100, it is NOT a leap year, unless;
% The year is also evenly divisible by 400. Then it is a leap year.


if ((mod(year,4) == 0) && ~(mod(year,100) == 0)) || (mod(year,400) == 0)
    nDays = 366;
else
    nDays = 365;
end

end