% Sctipt to run a simple PPP algorithm using the xona-pnt repo

%% inputs
% RINEX v3 observation file
filenameGnss = 'D:\PNT Data\Roof logs\swift-gnss-20200312-093212.sbp.obs';
% filenameGnss = 'C:\Users\kazgu\Desktop\Stanford\swift-gnss-20200312-093212.sbp.obs';

% need a configutation file to set where to put downloaded products.  The
% default included is called default.ini
configFile = 'config.ini';

% Truth position of the Xona roof antenna
truePosEcef = [-2706115.1823 -4278731.1983 3866392.5504];

% Three letter code indicating desired IGS analysis center
igsAc = 'GRG';

%%
% this should be multi-constellation!
constUse = [1 1 1 0 0];  % GPS | GLO | GAL | BDS | QZSS

% Initialize the filter
filter = navsu.estimators.pppFilter;
% filter = navsu.estimators.leastSq;

filter.PARAMS.states.RX_DCB_GLO = false;
filter.PARAMS.Q.POS = 0;
filter.PARAMS.Q.VEL = 0;

filter.PARAMS.measMask.f1 = [0 0 0]';
filter.PARAMS.measMask.f2 = [0 0 0]';
filter.PARAMS.measMask.f3 = [0 0 0]';

%% Read observation file
if ~exist('obsStruc','var')   
    disp('Reading observation file')
    [obsStruc, constellations, epochs, date, pos, interval, antoff, antmod,...
        rxmod] = navsu.readfiles.loadRinexObs(filenameGnss,'constellations',...
        navsu.readfiles.initConstellation(constUse(1),constUse(2),constUse(3),constUse(4),constUse(5)));
    
    obsGnssRaw.meas      = obsStruc;
    obsGnssRaw.PRN       = constellations.PRN;
    obsGnssRaw.constInds = constellations.constInds;
    obsGnssRaw.epochs    = epochs;
    obsGnssRaw.tLock     = [];
end

epochStart = min(epochs)+200*60;
downsampleFac = 30;

%%
if ~exist('corrData','var')
    disp('Loading corrections')
    % looking for most recent products
    jdRange0 = floor([navsu.time.epochs2jd(min(epochs)) navsu.time.epochs2jd(max(epochs))]+0.5)-0.5;
    jdRangeProd = (min(jdRange0)-1):(max(jdRange0)+1);
    
    [doyProd,yearProd] = navsu.time.jd2doy(jdRangeProd);
    doyProd = floor(doyProd);
    
    corrData = navsu.svOrbitClock('configFile',configFile,'constUse',constUse);
    
    corrData.settings.gpsEphCenter = igsAc;
    corrData.settings.gpsClkCenter = igsAc;
    corrData.settings.gloEphCenter = igsAc;
    corrData.settings.gloClkCenter = igsAc;
    corrData.settings.galEphCenter = igsAc;
    corrData.settings.galClkCenter = igsAc;
    
    % load MGEX ephemeris data into our correction object
    corrData.initOrbitData(yearProd,doyProd);
    
    % load MGEX clock data into our correction object
    corrData.initClockData(yearProd,doyProd);
    
    % load antenna phase center data for satellites
    filenameAtx = 'igs14_sats_only.atx';
    
    corrData.initAtxData(filenameAtx);
    
    % dcb products
    corrData.settings.dcbSource = 6;
    corrData.initDcb(yearProd,doyProd);
    
    % ionospheric data
    corrData.initIonoData(yearProd,doyProd);
end

%% preprocess observations
obsInds    = repmat([1 2 3 4],1,3); % 1 = code, 2 = carrier, 3 = snr, 4 = doppler
signalInds = kron(1:3,ones(1,4));
%           1                               2                               3                               4                                5
obsDes  = {{'C1C'} {'L1C'} {'S1C'}  {'D1C'} {'C2S'} {'L2S'} {'S2S'} {'D2S'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
    {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
    {'C1B'} {'L1B'} {'S1B'} {'D1B'} {'C7I'} {'L7I'} {'S7I'} {'D7I'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
    {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
    {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}  };


% signal pairs to include as iono-free combinations
ifPairs = [1 3;
    1 2];

[obsGnssi, dcbCorr0] = navsu.ppp.preprocessGnssObs(obsGnssRaw,obsInds,signalInds,...
    obsDes,ifPairs,corrData,'downsampleFac',downsampleFac,'epochStart',epochStart);

%% do the ppp lol
outData = navsu.ppp.runPpp(filter,obsGnssi,corrData);

%%
close all;

%% plot the skyplot
% figure;
% navsu.geo.skyPlot(outStruc.gnssData.az'*180/pi,outStruc.gnssData.el'*180/pi);

%%
filter.plotOutput(outData,'truePosEcef',truePosEcef);





