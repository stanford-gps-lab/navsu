% Sctipt to run a simple PPP algorithm using the navsu repo

clear variables;
close all;

%% inputs
% RINEX v3 observation file
% path to the obs file, can be copied to desired location from the
% navsu/examples folder.
filenameGnss = 'swift-gnss-20200312-093212.obs';
% Truth position of the Xona roof antenna
% (matches swift-gnss-20200312-093212.obs)
truePosEcef = [-2706115.1823 -4278731.1983 3866392.5504];

% run with precise or broadcast orbits?
usePreciseOrbits = false;

% use carrier smoothing?
useCarrierSmoothing = true;

% use the Klobuchar model to estimate the iono delay? (the use of dual
% frequency, iono-free measurements is unaffected by this)
useIonoModel = true;

% use the UNB3 model to estimate the tropospheric delay?
useTropoModel = true;

% set an elevation mask angle (default is 15 deg)
maskAngle = 5*pi/180; % in radians


% select frequencies to use
% two frequencies can be used for iono-free combination.
% Will be used in the order L1, L2, L5
freqs = [1 2]; % alternative = 1; = 2;

% which constellations should be used for positioning?
useConst = [1 1 1 0 0]; % [GPS GLO GAL BDS QZSS]

%% Read Rinex observation file

disp('Reading observation file')
[obsStruc, constellations, epochs] = navsu.readfiles.loadRinexObs(filenameGnss);
obsGnssRaw.meas      = obsStruc;
obsGnssRaw.PRN       = constellations.PRN;
obsGnssRaw.constInds = constellations.constInds;
obsGnssRaw.epochs    = epochs;
obsGnssRaw.tLock     = [];
% NOTE the absence of a lock time will inhibit carrier smoothing :-(

%% Read orbit data
% At least one file location must be specified to set up the orbit
% products. It is recommended to create a local copy of this example script
% and modify it accordingly.

% A configuration file should be specified- this tells the software where
% to place all of the downloaded products.  It can be modeled after the
% default.ini file that already exists. 
configFile = 'config.ini';

% Username and password file for NASA data/products download. THESE ARE NO
% LONGER REQUIRED, but their use can result in faster downloads.
% For their use see: 
% [1] https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
% [2] https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html
netrcFile = 'C:\...\cddislogin.netrc'; % change to your directory
% cddis also requires a cookie file when using cURL, see [1]
cookieFile = 'C:\...\nasa_cddis_cookies.txt'; % change to your directory

% Three letter code indicating desired IGS analysis center
igsAc = 'GRG';

disp('Loading corrections')
% looking for most recent products
jdRange0 = floor([navsu.time.epochs2jd(min(epochs)) ...
                  navsu.time.epochs2jd(max(epochs))]+0.5) - 0.5;
jdRangeProd = (min(jdRange0)-1):(max(jdRange0)+1);

[doyProd,yearProd] = navsu.time.jd2doy(jdRangeProd);
doyProd = floor(doyProd);

corrData = navsu.svOrbitClock('configFile', configFile, ...
                              'constUse', useConst, ... % defaults to GPS only
                              'netrcFile', netrcFile, ... % not required
                              'cookieFile', cookieFile); % not required
% can also be called as:
% corrData = navsu.svOrbitClock('configFile', configFile, ...
%                               'constUse', useConst);

corrData.settings.gpsEphCenter = igsAc;
corrData.settings.gpsClkCenter = igsAc;
corrData.settings.gloEphCenter = igsAc;
corrData.settings.gloClkCenter = igsAc;
corrData.settings.galEphCenter = igsAc;
corrData.settings.galClkCenter = igsAc;

% load MGEX precise ephemeris data into our correction object
% Only necessary if nav engine is to be run with PRECISE orbits
if usePreciseOrbits
    corrData.initOrbitData(yearProd, doyProd);

    % load MGEX clock data into our correction object
    % Only necessary if nav engine is to be run with PRECISE orbits
    corrData.initClockData(yearProd, doyProd);

    % load antenna phase center data for satellites
    % Only necessary if nav engine is to be run with PRECISE orbits
    filenameAtx = 'igs14_sats_only.atx';
    corrData.initAtxData(filenameAtx);
else
    % Load broadcast data
    % Only necessary if nav engine is to be run with BROADCAST orbits
    corrData.initBroadcastData(yearProd, doyProd);
end

%% Set up and run nav engine
% the navigation engine recognizes automatically if precise orbit data is
% provided in corrData. Prioritizes using precise data if it is available.
navEngine = navsu.lsNav.DFMCnavigationEngine(corrData);

% Alternatively can be called like this:
% navEngine = navsu.lsNav.DFMCnavigationEngine;
% navEngine.satEph = corrData;

% tell the engine whether to attempt carrier smoothing
navEngine.useCarrierSmoothing = useCarrierSmoothing;

% tell the engine whether to model the single frequency iono delay
navEngine.useIonoModel = useIonoModel;

% tell the engine whether to model the tropospheric delay
navEngine.useTropoModel = useTropoModel;

% OPTIONAL: Change elevation mask angle (stored in rad)
navEngine.elevMask = maskAngle; % default is 15 deg

% Due to the lack of a carrier lock time, no carrier smoothing will be
% performed on this dataset. If it were provided, the following parameters
% could be changed:

% OPTIONAL: change smoothing constants
navEngine.smoothingConstant = 300; % single frequency default is 100 sec
navEngine.smoothingConstantIF = 3000; % iono-free default is 1800 sec

% OPTIONAL: turn off carrier smoothing altogether
navEngine.useCarrierSmoothing = false; % default is true

% the below commented lines are summarized in the following call:
[posECEF, velECEF, tBias, R, prr, P, chi2stat, dop, useDF] = ...
    navEngine.batchPvtSolution(obsGnssRaw, useConst, freqs);

% select the satellites to use
% sats = ismember(obsGnssRaw.constInds, find(useConst));
% 
% % initialize all the outputs
% nEpochs = length(obsGnssRaw.epochs);
% posECEF = NaN(3, nEpochs);
% velECEF = NaN(3, nEpochs);
% tBias = NaN(sum(useConst), nEpochs);
% prr = NaN(sum(sats), nEpochs); % pseudorange residuals
% R = NaN(3+sum(useConst), 3+sum(useConst), nEpochs); % cov. matrix
% P = NaN(sum(sats), sum(sats), nEpochs); % residuals information matrix
% chi2stat = NaN(1, nEpochs); % chi-squared statistic
% 
% for ep = 1:length(obsGnssRaw.epochs)
%     [obsData, satIds] = navEngine.readRinexData(obsGnssRaw, sats, ep);
%     % limit to selected frequencies
%     obsData = structfun(@(x) x(:, freqs), obsData, 'UniformOutput', false);
%     
%     % do PVT computation
%     [posECEF(:, ep), tBias(:, ep), R(:, :, ep), prr(:, ep), P(:, :, ep)] = ...
%         navEngine.positionSolution(satIds, obsGnssRaw.epochs(ep), obsData);
%     
%     velECEF(:, ep) = ...
%         navEngine.velocitySolution(satIds, obsGnssRaw.epochs(ep), obsData);
%     
%     % also compute chi square statistic
%     s = isfinite(prr(:, ep));
%     chi2stat(ep) = prr(s, ep)' * P(s, s, ep) * prr(s, ep);
%                                             
% end

%% Plot position error
xyz = 'xyz';
figure;
for ii = 1:3
    subplot(3, 1, ii); hold on; grid on;
    % plot error
    plot(obsGnssRaw.epochs - obsGnssRaw.epochs(1), ...
         posECEF(ii, :) - truePosEcef(ii), ...
        'LineWidth', 2);
    % plot 2 sigma
    plot(obsGnssRaw.epochs - obsGnssRaw.epochs(1), ...
         2 * [1, -1] .* sqrt(squeeze(R(ii, ii, :))), ...
        'k--');
    ylabel(['ECEF ', xyz(ii), ' (m)'], ...
        'FontSize', 16, 'Interpreter', 'latex');
end
xlabel('Time (sec)', ...
       'FontSize', 16, 'Interpreter', 'latex');
sgtitle('ECEF position error and $2\sigma$ bounds (m)', ...
        'FontSize', 16, 'Interpreter', 'latex', 'FontWeight', 'bold');