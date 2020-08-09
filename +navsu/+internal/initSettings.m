function settings = initSettings(varargin)

% Parse inputs
p = inputParser;
p.addParameter('configFile', []); 
parse(p, varargin{:});
res = p.Results;
iniFile = res.configFile;

%% Pull info from the .ini file if available
if ~isempty(iniFile)
    iniData = navsu.thirdparty.ini2struct(iniFile);
    preciseProdDir = iniData.preciseproddir;
    obsDir         = iniData.obsdir;
else
    cdi = cd;
    cdi(strfind(cdi,'\')) = '/';
    preciseProdDir = [cdi '/data/'];
    obsDir         = [cdi '/data/'];
end


%% Data directories
% GPS----------------------------------------------------------------------
% Daily precise products (orbit and clock)
settings.preciseProdDir = [preciseProdDir 'precise-daily/'];
% Observations from IGS stations
settings.mgxObsDir      = [obsDir 'mgex-obs/'];
% High rate observations
settings.mgxHrObsDir    = [obsDir 'mgex-hr-obs/'];
% Temporary directory
settings.tempDir        = [preciseProdDir 'temp/'];
% Differential code biases
settings.dcbDir         = [preciseProdDir 'dcb/'];
% GPS navigation data from IGS stations
settings.navGpsDir      = [preciseProdDir 'nav-daily/gps/'];
% Multi-GNSS navigation data from IGS stations
settings.navMgxDir         = [preciseProdDir 'nav-daily/mgex/'];


%% Precise clock and orbit centers
% Precise ephemeris source
settings.gpsEphCenter   = 'com'; % NGA/IGS (final)/com/grm/wum
settings.gloEphCenter   = 'com'; % iac/com/grm
settings.galEphCenter   = 'com'; % com/grm/wum
settings.bdsEphCenter   = 'com'; % com/grm/wum

settings.highRate       = 1;     % Flag indicating use of high rate clock products

% Precise clock source
settings.gpsClkCenter   = 'com'; % MGEX/IGS source of high rate precise clock
settings.gloClkCenter   = 'com'; % MGEX/IGS source of high rate precise clock
settings.galClkCenter   = 'com'; % MGEX/IGS source of high rate precise clock
settings.bdsClkCenter   = 'com'; % MGEX/IGS source of high rate precise clock

settings.navSource   = 'sugl'; % sugl, iac
settings.MGEXSource = 's'; % THIS IS ACTUALLY GAL NAV SOURCE

% old settings still needed :(
settings.GloPephSource = 'MGEX';
settings.GloPclkSource = 'MGEX';
settings.nPolyFit       = 12;    % number of nominal points to use for polynominal fit
settings.pfit           = 8;     % order of fit for position interpolation
settings.cdfit          = 2;     % order of fit for clock interpolation

%% IGS Analysis Center settings
settings.orbCenter = {settings.gpsEphCenter ,settings.gloEphCenter ,settings.galEphCenter ,settings.bdsEphCenter ,'com'};
settings.clkCenter = {settings.gpsClkCenter ,settings.gloClkCenter ,settings.galClkCenter ,settings.bdsClkCenter ,'com'};

settings.orbitInterpMethod = 'lagrange';

%% Primary directory
baseDir = preciseProdDir;
tempDir = [baseDir filesep 'Temp Dir' filesep];

settings.tempDir = tempDir;

end






















