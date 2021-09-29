% Sctipt to run a simple PPP algorithm using the navsu repo

clear variables;
close all;

%% inputs
% RINEX v3 observation file
% path to the obs file, can be copied to desired location from the
% navsu/examples folder.
filenameGnss = '...\swift-gnss-20200312-093212.obs';
% Truth position of the Xona roof antenna
% (matches swift-gnss-20200312-093212.obs)
truePosEcef = [-2706115.1823 -4278731.1983 3866392.5504];

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

% Three letter code indicating desired IGS analysis center
igsAc = 'GRG';

% Define constellations to use (multi-constellation capable)
constUse = [1 1 1 0 0];  % GPS | GLO | GAL | BDS | QZSS

%% Read Rinex observation file

disp('Reading observation file')
[obsStruc, constellations, epochs, date, pos, interval, antoff, antmod,...
    rxmod] = navsu.readfiles.loadRinexObs(filenameGnss,'constellations',...
    navsu.readfiles.initConstellation(constUse(1),constUse(2),constUse(3),constUse(4),constUse(5)));

obsGnssRaw.meas      = obsStruc;
obsGnssRaw.PRN       = constellations.PRN;
obsGnssRaw.constInds = constellations.constInds;
obsGnssRaw.epochs    = epochs;
obsGnssRaw.tLock     = [];

epochStart = min(epochs);
downsampleFac = 10;

%% Initialize orbit propagation and corrections object
disp('Loading corrections')
% looking for most recent products
jdRange0 = floor([navsu.time.epochs2jd(min(epochs)) ...
                  navsu.time.epochs2jd(max(epochs))]+0.5) - 0.5;
jdRangeProd = (min(jdRange0)-1):(max(jdRange0)+1);

[doyProd,yearProd] = navsu.time.jd2doy(jdRangeProd);
doyProd = floor(doyProd);

corrData = navsu.svOrbitClock('configFile', configFile, ...
                              'constUse', constUse, ...
                              'netrcFile', netrcFile, ...
                              'cookieFile', cookieFile);

corrData.settings.gpsEphCenter = igsAc;
corrData.settings.gpsClkCenter = igsAc;
corrData.settings.gloEphCenter = igsAc;
corrData.settings.gloClkCenter = igsAc;
corrData.settings.galEphCenter = igsAc;
corrData.settings.galClkCenter = igsAc;

% load MGEX precise ephemeris data into our correction object
corrData.initOrbitData(yearProd,doyProd);

% load MGEX clock data into our correction object
corrData.initClockData(yearProd,doyProd);

% Load broadcast data
corrData.initBroadcastData(yearProd,doyProd);

% load antenna phase center data for satellites
filenameAtx = 'igs14_sats_only.atx';
corrData.initAtxData(filenameAtx);

% dcb products
corrData.initDcb(yearProd,doyProd);

% ionospheric data
corrData.initIonoData(yearProd,doyProd);

%% preprocess observations
[gnssMeas, dcbCorr0] = navsu.ppp.preprocessGnssObs(obsGnssRaw,...
     corrData,'downsampleFac',downsampleFac,'epochStart',epochStart,...
     'epochEnd',epochStart+30*60);
 
 % Build position measurement
% posMeas = navsu.ppp.buildPosMeas(gnssMeas.epochs,truePosEcef',0.1,1);

% velMeas = navsu.ppp.buildVelMeas(gnssMeas.epochs,[0 0 0]',0.01,1);


%% Estimate!
% Initialize the filter
filter = navsu.estimators.pppFilter;
% filter = navsu.estimators.leastSq;

filter.PARAMS.states.RX_DCB_GLO = false;
filter.PARAMS.Q.PXOS = 0;
filter.PARAMS.Q.VEL = 0;

filter.PARAMS.measMask.f1 = [0 0 1]';
filter.PARAMS.measMask.f2 = [0 0 0]';
filter.PARAMS.measMask.f3 = [0 0 0]';

filter.PARAMS.measUse.noVertVel = 0;


% choose precise or broadcast orbits (PRECISE is default)
% corrData.orbMode = 'BROADCAST';
% corrData.clkMode = 'BROADCAST';

corrData.orbMode = 'PRECISE';
corrData.clkMode = 'PRECISE';

outData = navsu.ppp.runPpp(filter,{gnssMeas },corrData);

%% Plot the results

filter.plotOutput(outData,'truthFile',truePosEcef);




