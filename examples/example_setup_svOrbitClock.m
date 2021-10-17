 
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
% NO LONGER NEEDED ON A MAC
netrcFile = 'C:\...\cddislogin.netrc'; % change to your directory
% cddis also requires a cookie file when using cURL, see [1]
cookieFile = 'C:\...\nasa_cddis_cookies.txt'; % change to your directory

% Initialize the GNSS correction data object
datai = navsu.svOrbitClock('constUse',   constUse, ...
                           'configFile', configFile, ... % defaults to GPS only
                           'netrcFile',  netrcFile, ... % not needed on Mac
                           'cookieFile', cookieFile); % not needed on Mac

% Year and day of year of interest- this is just arbitrary for this
% example
year = 2020;
doy  = 110;

%% Choose an example time frame and satellite
% Choose a time from the middle of the day (need to be able to interpolate
% using data from both sides)

% Build start time from day of year information 
epochStart = navsu.time.jd2epochs(navsu.time.doy2jd(year,doy)+0.5);  
dt = 10;  % Time interval in seconds
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
% run_ppp_multiconst.m example.)
datai.initDcb(year, doy);

% Initialize antenna phase center data
% (Accounts for offset between satellite COM and actual antenna phase
% center when propagating precise orbits.)
datai.initAtxData('igs14_sats_only.atx');

% Also initialize broadcast orbit data
% (In case we want to compute satellite positions using the broadcast nav
% files.)
datai.initBroadcastData(year, doy);


%% Interpolate the precise orbit and clock data for one satellite
disp(['Interpolating orbit and clock corrections for PRN ', ...
      num2str(prn), ' to the desired epochs'])
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
[svPosB, svVelB] = datai.propagate(prnVec, constIndVec, epochs);
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


%% Interpolate the precise orbit and clock data for ALL satellite
disp('Interpolating orbit and clock corrections for all satellites.')
% Orbit data
% This returns antenna phase center positions since we initialized the ATX
% data
datai.orbMode = 'PRECISE'; % should be the default
datai.clkMode = 'PRECISE'; % should be the default
% preallocate
[svPos, svVel, svPosB, svVelB] = deal(NaN(length(datai.PEph.PRN), 3, length(epochs)));
[svBias, svBiasB] = deal(NaN(length(datai.PEph.PRN), length(epochs)));

for ep = 1:length(epochs)
    % pos, vel
    [svPos(:, :, ep), svVel(:, :, ep)] = datai.propagate(datai.PEph.PRN, ...
                                                   datai.PEph.constellation, ...
                                                   epochs(ep)*ones(size(datai.PEph.PRN)));
    % Clock data
    svBias(:, ep) = datai.clock(datai.PEph.PRN, ...
                             datai.PEph.constellation, ...
                             epochs(ep)*ones(size(datai.PEph.PRN)));
end

% Now propagate the orbits using the broadcast data as a comparison
datai.orbMode = 'BROADCAST'; % instead of the default 'PRECISE'
datai.clkMode = 'BROADCAST'; % instead of the default 'PRECISE'

for ep = 1:length(epochs)
    [svPosB(:, :, ep), svVelB(:, :, ep)] = datai.propagate(datai.PEph.PRN, ...
                                                     datai.PEph.constellation, ...
                                                     epochs(ep)*ones(size(datai.PEph.PRN)));
    svBiasB(:, ep) = datai.clock(datai.PEph.PRN, ...
                                 datai.PEph.constellation, ...
                                 epochs(ep)*ones(size(datai.PEph.PRN)));
end

svPosErr = squeeze(sqrt(sum((svPosB - svPos).^2, 2)));
svVelErr = squeeze(sqrt(sum((svVelB - svVel).^2, 2)));
svClkErr = (svBiasB - svBias) * navsu.constants.c;

%% Analyze the error of the broadcast orbits
constNames = {'GPS'; 'GLO'; 'GAL'; 'BDS'; 'QZSS'; 'SBAS'};
for consti = find(unique(datai.PEph.constellation))'
    % plot orbit position error
    figure; hold on; grid on;
    plot(epochs - epochs(1), ...
         svPosErr(datai.PEph.constellation == consti, :), ...
         'LineWidth', 2);
    xlabel('Time (sec)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel([constNames{consti}, ' absolute orbit error (m)'], ...
        'FontSize', 18, 'Interpreter', 'latex');
    
    % plot velocity error
    figure; hold on; grid on;
    plot(epochs - epochs(1), ...
         svVelErr(datai.PEph.constellation == consti, :), ...
         'LineWidth', 2);
    xlabel('Time (sec)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel([constNames{consti}, ' absolute velocity error (m/s)'], ...
        'FontSize', 18, 'Interpreter', 'latex');
    
    % and clock error
    figure; hold on; grid on;
    plot(epochs - epochs(1), ...
         svClkErr(datai.PEph.constellation == consti, :), ...
         'LineWidth', 2);
    xlabel('Time (sec)', 'FontSize', 18, 'Interpreter', 'latex');
    ylabel([constNames{consti}, ' clock error (m)'], ...
        'FontSize', 18, 'Interpreter', 'latex');
end
