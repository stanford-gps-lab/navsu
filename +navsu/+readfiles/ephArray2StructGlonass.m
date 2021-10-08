function eph = ephArray2StructGlonass(array, filename, leapSecond)
% ephArray2StructGlonass
% DESCRIPTION:
%   Converts navigation message arrays from GLONASS to a structure. 
% INPUT:
%   array      - parsed navigation message data from
%                navsu.readfiles.loadRinexNav
%   filename   - name of the RINEX file that the data came from
%   leapSecond - number of leap seconds read from the header of the RINEX
%                file
%
% OUTPUT:
%   eph        - structure containing all of the navigation message data
%
% See also: navsu.readfiles.loadRinexNav

if size(array, 2) < 22
    eph = [];
    return
end

y = array(:, 2);
z = (y >= 80) & (y <= 99);
array(z, 2) = array(z, 2) + 1900;
z = (y >= 0) & (y <= 79);
array(z, 2) = array(z, 2) + 2000;
if any(y > 99) || any(y < 0)
%     fprintf(2, 'Warning: Perhaps incorrect years range from %d to %d', min(y), max(y));
end

% eph.filename         = filename;
% eph.leapSecond       = leapSecond;
% eph.PRN              = array(:, 1);
% [eph.GPS_week_num, eph.Toc, eph.GPS_weekday] = navsu.time.utc2gps(array(:, 2:7), false);
% eph.clock_bias       = array(:, 8);
% eph.clock_drift      = array(:, 9);
% eph.Tframe           = array(:, 10);
% eph.X                = array(:, 11);
% eph.Xd               = array(:, 12);
% eph.Xdd              = array(:, 13);
% eph.health           = array(:, 14);
% eph.Y                = array(:, 15);
% eph.Yd               = array(:, 16);
% eph.Ydd              = array(:, 17);
% eph.freqN            = array(:, 18);
% eph.Z                = array(:, 19);
% eph.Zd               = array(:, 20);
% eph.Zdd              = array(:, 21);
% eph.AOI              = array(:, 22);

ToE = datenum(array(:, 2:7));
% nd = mode(floor(ToE)) - datenum([1980 1 6]);
nd2 = floor(ToE)- datenum([1980 1 6]);
eph.GPS_week_num     = floor(nd2 / 7);
eph.GPS_weekday      = mod(nd2, 7);
eph.ToE              = round((ToE - datenum([1980 1 6])) * 86400);   % seconds since 1980-1-6
eph.filename         = filename;
eph.leapSecond       = leapSecond;
eph.PRN              = array(:, 1);
eph.clock_bias       = array(:, 8);
eph.clock_drift      = array(:, 9);
eph.tk               = array(:, 10);
eph.P                = array(:, [11 15 19]) .* 1e3;
eph.V                = array(:, [12 16 20]) .* 1e3;
eph.A                = array(:, [13 17 21]) .* 1e3;
eph.health           = array(:, 14);
eph.freqNum          = array(:, 18);
eph.AoOper           = array(:, 22);