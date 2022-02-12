function eph = ephArray2Struct(array, filename, leapSecond, constellation)
% ephArray2Struct
% DESCRIPTION:
%   Converts navigation message arrays from GPS, Galileo, or Beidou to a
%   structure.
% INPUT:
%   array      - parsed navigation message data from
%                navsu.readfiles.loadRinexNav
%   filename   - name of the RINEX file that the data came from
%   leapSecond - number of leap seconds read from the header of the RINEX
%                file
%   constellation - name of the constellation associated with the data here
%                'GPS','GAL','BDS'
% OUTPUT:
%   eph        - structure containing all of the navigation message data
%
% See also: navsu.readfiles.loadRinexNav

if (size(array, 2) < 36) && (size(array, 2) ~= 22) % second term is for SBAS
%     error('Incorrect size of ephemerides array: %d', size(array, 2));
    eph = [];
    return
end

if nargin < 4
    constellation = 'GPS';
end

y = array(:, 2);
z = (y >= 80) & (y <= 99);
array(z, 2) = array(z, 2) + 1900;
z = (y >= 0) & (y <= 79);
array(z, 2) = array(z, 2) + 2000;
% if any(y > 99) || any(y < 0)
%     fprintf(2, 'Warning: Perhaps incorrect years range from %d to %d', ...
%             min(y), max(y));
% end

eph.filename         = filename;
eph.leapSecond       = leapSecond;
eph.PRN              = array(:, 1);
[~, eph.Toc, eph.GPS_weekday] = navsu.time.jd2gps( ...
    navsu.time.cal2jd(array(:,2), array(:,3), ...
        array(:,4)+array(:,5)/24+array(:,6)/(24*60)+array(:,7)/(86400)));
eph.clock_bias       = array(:, 8);

if strcmp(constellation,'SBAS')
    % rows 11-22 mostly follow GLONASS(!) record layout...
    eph.frequency_bias = array(:, 9);
    eph.TTOM           = array(:, 10);
    eph.X              = array(:, 11);
    eph.Xd             = array(:, 12);
    eph.Xdd            = array(:, 13);
    eph.health         = array(:, 14);
    eph.Y              = array(:, 15);
    eph.Yd             = array(:, 16);
    eph.Ydd            = array(:, 17);
    eph.URA            = array(:, 18);
    eph.Z              = array(:, 19);
    eph.Zd             = array(:, 20);
    eph.Zdd            = array(:, 21);
    eph.IODN           = array(:, 22);
    return
end

eph.clock_drift      = array(:, 9);
eph.clock_drift_rate = array(:, 10);
% if strcmp(constellation,'GAL')
%     % IODE is not defined for GAL -- the most closely analogous quantity is
%     % called IODNav -- so take this value from ToE (see below) instead.
%     %
%     % QUESTION: per the (somewhat hard-to-follow) discusslion presented in
%     % <https://destevez.net/2019/09/ephemeris-quality-during-the-galileo-outage/>,
%     % might it be more appropriate to take this value from eph.TTOM (that
%     % is, array(:, 35)) instead?
%     eph.IODE         = array(:, 19); % this is just eph.Toe
% else
eph.IODE             = array(:, 11); % IODN in the case of Galileo
% end
eph.Crs              = array(:, 12);
eph.Delta_n          = array(:, 13);
eph.M0               = array(:, 14);
eph.Cuc              = array(:, 15);
eph.e                = array(:, 16);
eph.Cus              = array(:, 17);
eph.sqrtA            = array(:, 18);
eph.Toe              = array(:, 19);
eph.Cic              = array(:, 20);
eph.OMEGA            = array(:, 21);
eph.Cis              = array(:, 22);
eph.i0               = array(:, 23);
eph.Crc              = array(:, 24);
eph.omega            = array(:, 25);
eph.OMEGA_DOT        = array(:, 26);
eph.I_DOT            = array(:, 27);
eph.codes_on_L2      = array(:, 28);
eph.GPS_week_num     = array(:, 29);
eph.L2_P_data_flag   = array(:, 30);

% Developer note- URA and SISA translation not yet available- ask Kaz :)
if strcmp(constellation, 'GPS')
%     eph.accuracy         = TranslateURA(array(:, 31));
elseif strcmp(constellation,'GAL')
%     eph.accuracy         = translateSISA(array(:, 31));
elseif strcmp(constellation,'BDS')
   eph.accuracy = array(:,31);
%    eph.aodClock = translateBdsAod(array(:,36));
%    eph.aodEph   = translateBdsAod(array(:,11));
end
eph.accuracy = array(:,31);

eph.health           = array(:, 32);
eph.TGD              = array(:, 33);
if strcmp(constellation,'GPS')
    eph.IODC         = array(:, 34);
else % no iodc saved? idk
    eph.IODC = eph.IODE;
end

if any(strcmp(constellation, {'GAL', 'BDS'}))
    % these have different group delays for different frequency
    % combinations
    eph.TGD2 = array(:, 34);
end

eph.TTOM             = array(:, 35);
eph.Fit_interval     = array(:, 36);

if size(array,2) == 38
   eph.suglInfo1 = array(:,37);
   eph.suglInfo2 = array(:,38);

end


