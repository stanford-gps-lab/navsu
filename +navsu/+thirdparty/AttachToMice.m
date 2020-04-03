function AttachToMice(miceDir)

% Initialize MICE utility for solar/planetary ephemeris and earth rotation
% parameters
if exist('src\mice\','dir')
   % don't need to attach- get out of here!
   return
end

addpath([miceDir '\src\mice\'])
addpath([miceDir '\lib\'])

% Planetary ephemerides
cspice_furnsh([miceDir '\data\de430.bsp'])

% Earth rotation
cspice_furnsh([miceDir '\data\earth_fixed.tf'])
cspice_furnsh([miceDir '\data\earth_assoc_itrf93.tf'])
cspice_furnsh([miceDir '\data\earth_latest_high_prec.bpc'])
cspice_furnsh([miceDir '\data\earth_070425_370426_predict.bpc'])
cspice_furnsh([miceDir '\data\earth_000101_191225_191003.bpc'])
cspice_furnsh([miceDir '\data\earth_000101_161009_160718.bpc'])

% Leap seconds
cspice_furnsh([miceDir '\data\naif0011.tls'])
cspice_furnsh([miceDir '\data\naif0012.tls'])

end