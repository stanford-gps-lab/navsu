 
% Input vector for navsu.svOrbitClock telling it what constellations to
% load.  First element is GPS, then GLONASS, Galileo, BDS, SBAS
% respectively.  
constUse = [1 1 0 0 0];

% A configuration file should be specified- this tells the software where
% to place all of the downloaded products.  It can be modeled after the
% default.ini file that already exists. 
configFile = 'config.ini';

% Username and password file for NASA data/products download. See: 
% [1] https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
% [2] https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html
netrcFile = 'C:\...\cddislogin.netrc'; % change to your directory
% cddis also requires a cookie file when using cURL, see [1]
cookieFile = 'C:\...\nasa_cddis_cookies.txt'; % change to your directory

% Initialize the GNSS correction data object
datai = navsu.svOrbitClock('constUse',   constUse, ...
                           'configFile', configFile, ...
                           'netrcFile',  netrcFile, ...
                           'cookieFile', cookieFile);

% Year and day of year of interest- this is just arbitrary for this
% example
year = 2020;
doy  = 110;

%% Choose an example time frame and satellite
% Choose a time from the middle of the day (need to be able to interpolate
% using data from both sides)

% Build start time from day of year information 
epochStart = navsu.time.jd2epochs(navsu.time.doy2jd(year,doy)+0.5);  
dt = 60;  % Time interval in seconds
epochs = (epochStart:dt:(epochStart+3600))'; % just running for one hour

% Choosing an arbitrary GPS or GLONASS satellite 
prn = 10;
constInd = 1;  % GPS = 1, GLONASS = 2, GAL = 3, ...

% the PRN, constInd, and epoch inputs must all be the same length to the
% orbit interpolation function.
% Could be run for multiple PRNs at once, too!
prnVec = prn*ones(size(epochs));
constIndVec = constInd*ones(size(epochs));

%% Initialize Orbit and clock corrections
disp('Downloading and/or parsing orbit and clock corrections')
% Download and/or parse the precise orbital data
datai.initOrbitData(year, doy);

% Download and/or parse the clock data
datai.initClockData(year, doy);

% Initialize ionospheric data
% (corrections used for single frequency positioning)
datai.initIonoData(year, doy);

% Initialize DCB data
% (corrections used only when processing actual measurements. See
% run_ppp_multiconst.m example.
datai.initDcb(year, doy);

% Initialize antenna phase center data
% (Accounts for offset between satellite COM and actual antenna phase
% center when propagating precise orbits.)
datai.initAtxData('igs14_sats_only.atx');

% Also initialize broadcast orbit data
% (In case we want to compute satellite positions using the broadcast nav
% files.)
datai.initBroadcastData(year, doy);


%% Interpolate the precise orbit and clock data
disp('Interpolating orbit and clock corrections to the desired epochs')
% Orbit data
% This returns antenna phase center positions since we initialized the ATX
% data
datai.orbMode = 'PRECISE'; % should be the default
[svPos,svVel] = datai.propagate(prnVec, constIndVec, epochs);

% Clock data
datai.clkMode = 'PRECISE'; % should be the default
svBias = datai.clock(prnVec, constIndVec, epochs);

% Now propagate the orbits using the broadcast data as a comparison
datai.orbMode = 'BROADCAST'; % instead of the default 'PRECISE'
[svPosB, SVVELB] = datai.propagate(prnVec, constIndVec, epochs);
datai.clkMode = 'BROADCAST'; % instead of the default 'PRECISE'
svBiasB = datai.clock(prnVec, constIndVec, epochs);

%% Analyze the error of the broadcast orbits
figure; hold on; grid on;
plot(epochs - epochs(1), ...
    [svPosB - svPos, (svBiasB - svBias)*navsu.constants.c], ...
    'LineWidth', 2);
xlabel('Time since start in (sec)', 'FontSize', 14, 'Interpreter', 'latex')
ylabel('ECEF Orbit error in (m)', 'FontSize', 14, 'Interpreter', 'latex')
legend('x', 'y', 'z', 'clock', ...
    'Location', 'best', 'FontSize', 14, 'Interpreter', 'latex')
