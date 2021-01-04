function settings = initSettings(varargin)

% Parse inputs
p = inputParser;
p.addParameter('configFile', []); 
p.addParameter('netrcFile',[]);
p.addParameter('cookieFile',[]);
p.addParameter('cookieDir',[]);
parse(p, varargin{:});
res = p.Results;
iniFile = res.configFile;

% cookieDir
%   Local folder where cookies from NASA CDDIS web access can/should be
%   stored.  
% cookieFile
%   A storage file can also be specified for cookie storage.
% netrcFile
%   Username and password file for NASA CDDIS data/product access.
%   Information about general access can be found at [1], and specifics of
%   what this file should look like can be found at [2] once a username and
%   password have been created. 

% [1] https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
% [2] https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html

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


%% NASA login/download setup
defaultCookieFilename = 'nasa_cddis_cookies.txt';
if ~isempty(res.cookieFile)
    % A specific filename of the cookie was provided
    settings.cookieFile = res.cookieFile; 
elseif ~isempty(res.cookieDir)
    % only a directory was provided- 
    settings.cookieFile = fullfile(res.cookieDir,defaultCookieFilename);
else
    % No location information about where to place teh cookie file was
    % provided
    settings.cookieFile = fullfile(preciseProdDir,defaultCookieFilename);
end
settings.netrcFile = res.netrcFile;

%% Data directories
settings.baseDir = preciseProdDir;
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
tempDir = [baseDir '\Temp Dir\'];

settings.tempDir = tempDir;

end






















