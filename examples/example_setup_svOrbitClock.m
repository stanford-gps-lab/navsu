 
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
netrcFile = 'C:\Users\kazgu\Documents\cddislogin.netrc';

% Initialize the GNSS correction data object
datai = navsu.svOrbitClock('constUse',[1 1 0 0 0],'netrcFile',netrcFile);

% Year and day of year of interest- this is just arbitrary for this
% example
year = 2020;
doy  = 180;

%% Choose an example time frame and satellite
% Choose a time from the middle of the day (need to be able to interpolate
% using data from both sides)

% Start time is built from day of year information 
epochStart = navsu.time.jd2epochs(navsu.time.doy2jd(year,doy)+0.5);  
dt = 60;  % Time interval in seconds
epochs = (epochStart:dt:(epochStart+3600))'; % just running for one hour

% Choosing an arbitrary GPS or GLONASS satellite 
prn = 10;
constInd = 1;  % GPS = 1, GLONASS = 2, GAL = 3

% the PRN, constInd, and epoch inputs must all be the same length to the
% orbit interpolation function
prnVec = prn*ones(size(epochs));
constIndVec = constInd*ones(size(epochs));

%% Initialize Orbit and clock corrections
disp('Downloading and/or parsing orbit and clock corrections')
% Download and/or parse the orbital data
datai.initOrbitData(year,doy);

% Download and/or parse the clock data
datai.initClockData(year,doy);

% Interpolate the precise orbit and clock data
disp('Interpolating orbit and clock corrections to the desired epochs')
% Orbit data
[svPos,svVel] = datai.propagate(prnVec,constIndVec,epochs);

% Clock data
svBias = datai.clock(prnVec,constIndVec,epochs);

%% Initialize ionospheric data
datai.initIonoData(year,doy);

%% Initialize DCB data
datai.initDcb(year,doy);

%% Initialize antenna phase center data
datai.initAtxData('igs14_sats_only.atx');























