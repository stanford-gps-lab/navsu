function [GPSweek, GPSsecond, GPSweekday, leapSeconds] = utc2gps(utc_time, useLeapSecond)% [GPSweek, GPSsecond] = utc2gps(utc_time, leapSecond)% Convert UTC times to the equivalent time in GPS weeks and GPS seconds%% Input:  %   utc_time - m-by-6 matrix containing m full date vectors respectively. A%              full date vector has six elements, specifying year, month, %              day, hour, minute, and second, in that order.%   useLeapSecond - boolean scalar. If "true" then leap seconds will be%                   counted in; if "false" or omitted then no leap seconds%                   will be counted in.% Output: %   GPSweek - m-by-1 vector%   GPSsecond - m-by-1 vector%% Functions called: None%% Notice: the leapSecondsTable in the codes should be updated if new leap%         second is introduced.% Written by: Liang Heng  08/01/2009% Copyright (c) 2009 by GPS Research Laboratory, Stanford UniversityleapSecondsTable = [    1981           7           1           0           0           0           1    1982           7           1           0           0           0           2    1983           7           1           0           0           0           3    1985           7           1           0           0           0           4    1988           1           1           0           0           0           5    1990           1           1           0           0           0           6    1991           1           1           0           0           0           7    1992           7           1           0           0           0           8    1993           7           1           0           0           0           9    1994           7           1           0           0           0          10    1996           1           1           0           0           0          11    1997           7           1           0           0           0          12    1999           1           1           0           0           0          13    2006           1           1           0           0           0          14    2009           1           1           0           0           0          15    2012           7           1           0           0           0          16    2015           7           1           0           0           0          17    2017           1           1           0           0           0          18];leapDayNum = datenum(leapSecondsTable(:, 1:6));utcDayNum = datenum(utc_time);len = length(utcDayNum);GPSweek = NaN(len, 1);GPSsecond = NaN(len, 1);GPSweekday = NaN(len, 1);leapSeconds = nan(len,1);if (nargin >= 2) && useLeapSecond    for i = 1 : len        ind = find(utcDayNum(i) >= leapDayNum, 1, 'last');        if isempty(ind)            fprintf(2, 'Warning: Unacceptable utc_time(%d, :) = ''%s''\n', i, datestr(utcDayNum(i), 0));            continue        end        SecondNum = etime(utc_time(i, :), [1980 1 6 0 0 0]) + leapSecondsTable(ind, 7);        GPSweek(i) = floor(SecondNum / 604800);        GPSsecond(i) = mod(SecondNum, 604800);        leapSeconds(i) = leapSecondsTable(ind,7);    endelse    for i = 1 : len        SecondNum = etime(utc_time(i, :), [1980 1 6 0 0 0]);        GPSweek(i) = floor(SecondNum / 604800);        GPSsecond(i) = mod(SecondNum, 604800);        GPSweekday(i) = floor(GPSsecond(i) / 86400);    endend