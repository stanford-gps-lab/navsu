function settings = initSettings(varargin)

% Parse inputs
p = inputParser;
p.addParameter('iniFile', 'config.ini'); 
parse(p, varargin{:});
res = p.Results;
iniFile = res.iniFile;

%% Pull info from the .ini file
iniData = navsu.thirdparty.ini2struct(iniFile);
miceDir        = iniData.micedir;
preciseProdDir = iniData.preciseproddir;
obsDir         = iniData.obsdir;

%% Initialize directories and other settings for broadcast ephemeris and
% clock analysis

%% Filter settings
% Minimum SNR for frequency 1
settings.PARAMS.snrMin1 = 0;
% Minimum SNR for frequency 22
settings.PARAMS.snrMin2 = 0;
% Elevation mask angle
settings.PARAMS.elMask  = 7.5; % 12.5 % degrees (obviously)

% Tropo model to use
settings.PARAMS.tropModel = 'GMF'; % 'GMF','GO','UNB3'

% Iono mode to use for single frequency
settings.PARAMS.ionoMode = 'TECSTATE'; % 'PARAM','TEC','LOSEST','TECSTATE'
% 'TECSTATE' uses a combination of a TEC map and dual frequency
% measurements to produce an estimate of the iono delay for each LOS

% number of subsets to use (number to exclude)
settings.PARAMS.nOutSubset = 1;
settings.PARAMS.nMaxSubset = 100;

settings.PARAMS.faultResetSubsets = false;

% reduced state simplified filter mode
settings.PARAMS.simplifiedMode = false;

% additional states
settings.PARAMS.states.mpCode = true;
settings.PARAMS.states.mpCarr = true;
settings.PARAMS.states.dcb    = true;
settings.PARAMS.states.brdc   = false;
settings.PARAMS.states.waas   = false;

% additional measurements
settings.PARAMS.meas.sameHeight = false;

% single frequency cycle slip check?
settings.PARAMS.cycleSlipSf = true;
settings.PARAMS.cycleSlipDet = 'GEOMETRYFREE';  % 'RECEIEVEROUTPUT'

% INS use settings
settings.PARAMS.INS = false;
settings.PARAMS.states.attitude = false;
settings.PARAMS.states.attitudeMode = 'EULER'; % 'QUATERNION','EULER'
settings.PARAMS.states.attitude_rate = false;

settings.PARAMS.states.acc_scale = false;
settings.PARAMS.states.w_scale   = false;

settings.PARAMS.reducedMeas = true;
settings.PARAMS.alternateMeasUpdate = false;

% Settings related to use of broadcast ephemeris
settings.PARAMS.orbMode = 'PRECISE'; % 'PRECISE' or 'BROADCAST'
settings.PARAMS.clkMode = 'PRECISE';
settings.PARAMS.iodControl = 0;
settings.PARAMS.useSbas    = 0;

% Number of satellites to group per subset
settings.PARAMS.subsetGroupSize = 1;

settings.PARAMS.R_ACC = 0.19*1e-3*sqrt(500);
settings.PARAMS.R_GYRO = 0.02*pi/180*sqrt(100);
settings.PARAMS.IMU_ARM = -[-0.262 0.944 1.575]';

% Initial covariance values
settings.PARAMS.SIGMA_RXB0      = 100^2;   % receiver clock
settings.PARAMS.SIGMA_RXBD0     = 1e4^2;   % receiver clock rate
settings.PARAMS.SIGMA_POS0      = 3^2;    % position
settings.PARAMS.SIGMA_VEL0      = 0.1^2;    % velocity
settings.PARAMS.SIGMA_ACC0      = 0.1^2;    % acceleration
settings.PARAMS.SIGMA_TROP0     = 0.1^2;   % troposphere
settings.PARAMS.SIGMA_AMB0      = 100^2;    % carrier phase ambiguities
settings.PARAMS.SIGMA_MP0       = 2^2;     % multipath (code, per LOS)
settings.PARAMS.SIGMA_MPPH0     = 0.05^2;  % multipath (carrier, per LOS)
settings.PARAMS.SIGMA_IONO0     = 10^2;    %  iono estimate given dual freq
settings.PARAMS.SIGMA_IONO0_SF  = 10^2;    % iono estimate given no only single freq
settings.PARAMS.SIGMA_RXDCB0    = 20^2;    % differential code bias
settings.PARAMS.SIGMA_DRXDCB0   = 1^2;     % frequency dependent differential code bias
settings.PARAMS.SIGMA_BRDC0     = 2.4^2;   % broadcast ephemeris error
settings.PARAMS.SIGMA_WAAS0     = 0.4^2;   % WAAS correction error

% INS related initial covariance
settings.PARAMS.SIGMA_ATTITUDE0 = (100*pi/180)^2;
settings.PARAMS.SIGMA_ACC_BIAS0 =(1e-0*70e-3)^2;
settings.PARAMS.SIGMA_W_BIAS0   = (1e-0*1*pi/180)^2;
settings.PARAMS.SIGMA_ATT_RATE0 = (1e-0*2*pi/180)^2;
settings.PARAMS.SIGMA_ACC_SCALE0 = (1e-0*0.04)^2;
settings.PARAMS.SIGMA_W_SCALE0   = (1e-0*0.01)^2;

% Process noise values
settings.PARAMS.Q_RXB       = 100^2;        % receiver clock
settings.PARAMS.Q_RXBD      = 1e4^2;        % receiver clock rate
settings.PARAMS.Q_POS       = 1^2;          % position
settings.PARAMS.Q_VEL       = 0.5^2;        % velocity
settings.PARAMS.Q_ACC       = 0.5^2;        % acceleration
settings.PARAMS.Q_TROP      = 0.002^2/3600; % troposphere
settings.PARAMS.Q_AMB_SF    = 0^2+0*0.02^2;          % carrier phase ambiguities
settings.PARAMS.Q_AMB_DF    = 0^2+0*0.02^2;          % carrier phase ambiguities
settings.PARAMS.Q_MP        = 0.25^2;       % multipath (code, per LOS, also scales with elevation)
settings.PARAMS.Q_MPPH      = (2*0.0008^2+0.02^2);   % multipath (carrier, per LOS)
settings.PARAMS.Q_IONO      = 0.005^2;      % ionosphere 0.005
settings.PARAMS.Q_IONOPARAMS = 0.0001^2;    % parameterized iono state
settings.PARAMS.Q_RXDCB     = 0^2;          % differential code bias
settings.PARAMS.Q_DRXDCB    = 0;            % frequency dependent differential code bias
settings.PARAMS.Q_BRDC      = 0.02^2;       % broadcast ephemeris error
 
settings.PARAMS.IONO_F0 = 1575.42e6;

settings.PARAMS.FLAG_Q_HEADING = 1;
settings.PARAMS.SIGMA_POS_FUR = [0.772   0.1077   1.186];
    
% Multipath state parameters
settings.PARAMS.ALPHA_MP0 = 0.97; %0.97
settings.PARAMS.TAU_MP0   = 100; %0.97
settings.PARAMS.TAU_MPPH0 = 150;

% Amount of time without measurements avaialable to elapse before state removeal
settings.PARAMS.tNoMeasRemove = 60*5;

switch settings.PARAMS.ionoMode
    case 'TEC'
        prNoiseBase = 100;
        phNoiseBase = 0.05;
    case 'PARAM'
        prNoiseBase = 10;
        phNoiseBase = 5;
    case 'LOSEST'
        prNoiseBase = 3;
        phNoiseBase = 5;
    case 'TECSTATE'
        prNoiseBase = 3;
        phNoiseBase = 0.03;
end

% Std of code phase measurements
settings.PARAMS.MEAS_PR    = [prNoiseBase*[1 2 1 1 1]; % single frequency
                              4 7 1 2 2]; % dual frequency
% Std of carrier phase measurements
settings.PARAMS.MEAS_CP    = [phNoiseBase*[1 1 1 1 1];
                              0.04 0.04 0.03 0.03 0.03 ]; %0.03
% Std of doppler measurements
settings.PARAMS.MEAS_DOP = [0.3 0.3 0.3 0.03 0.03];
% settings.PARAMS.MEAS_DOP = [0.03 0.03 0.03 0.03 0.03];

% Cycle slip threshold for time differenced geometry-free combination
settings.PARAMS.slipThresh = 0.15;% 0.05; %20

% Measurement exclusion threshold for using simple residual check
settings.PARAMS.prExcludeResidThresh = 10;
settings.PARAMS.phExcludeResidThresh = 0.1;
settings.PARAMS.dopExcludeResidThresh = 10*settings.PARAMS.MEAS_DOP; % 3.89x
settings.PARAMS.posSeedExcludeResidThresh = 100;
% TEC nComp measurement noise
settings.PARAMS.TECmeasStd = 0.5; % [m] 0.05

% Satellite fault rate
settings.PARAMS.Psat = 1e-5;

% "noise" connecting single frequency to dual frequency parameters
settings.PARAMS.singleToDualFreqStd = 0.01;

% subset position closeness threshold
settings.PARAMS.subsetDistanceExclude = 0.001;

% Speed of light
settings.PARAMS.c =299792458;

% Desired observations (Rows are for GPS, GLONASS, Galileo, BDS
% (unsupported) and QZSS (unsupported)
%           F1 CODE             F1 CARRIER         F1 SNR              F2 CODE               F2 CARRIER          F2 SNR              F1 DOPPLER          F2 DOPPLER
obsDes  = {{'C1C'}             {'L1C'}             {'S1C'}             {'C2P' 'C2W' 'C2X'}  {'L2P' 'L2W' 'L2X'} {'S2P' 'S2W' 'S2X'} {'D1C'}             {'D2P' 'D2W' 'D2Y'}% G
           {'C1C' 'C1C'}       {'L1C' 'L1C'}       {'S1C' 'S1C'}       {'C2C' 'C2W' 'C2C'}  {'L2C' 'L2W' 'L2C'} {'S2C' 'S2W' 'S2C'} {'D1C' 'D1X' 'D1C'} {'D2C' 'D2X' 'D2C'}% R
           {'C1B' 'C1X' 'C1C'} {'L1B' 'L1X' 'L1C'} {'S1B' 'S1X' 'S1C'} {'C8X' 'C5Q' 'C5X'}  {'L8X' 'L5Q' 'L5X'} {'S8X' 'S5Q' 'S5X'} {'D1B' 'D1X' 'D1C'} {'D5I' 'D5Q' 'D5X'}% E
           {'C1C'}             {'L1C'}             {'S1C'}             {'C2P' 'C2W' 'C2Y'}  {'L2P' 'L2W' 'L2Y'} {'S2P' 'S2W' 'S2Y'} {'D1B' 'D1X' 'D1C'} {'D2B' 'D2X' 'D2C'}% C
           {'C1C'}             {'L1C'}             {'S1C'}             {'C2P' 'C2W' 'C2Y'}  {'L2P' 'L2W' 'L2Y'} {'S2P' 'S2W' 'S2Y'} {'D1B' 'D1X' 'D1C'} {'D2B' 'D2X' 'D2C'}}; % S

settings.obsDes = obsDes;

%% Protection level settings
settings.PARAMS.pl_ss_exact_optfa 	= 0; 
settings.PARAMS.pl_ss_approx1_optfa = 1;
settings.PARAMS.pl_ss_approx2_optfa = 0;
settings.PARAMS.pl_ss_exact 		= 0;
settings.PARAMS.pl_ss_approx1 		= 0;
settings.PARAMS.pl_ss_approx2 		= 0;
settings.PARAMS.pl_chi2_exact 		= 0;
settings.PARAMS.pl_chi2_approx1 	= 0; 
settings.PARAMS.pl_chi2_approx2 	= 0;


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

%% Filter settings
settings.orbCenter = {settings.gpsEphCenter ,settings.gloEphCenter ,settings.galEphCenter ,settings.bdsEphCenter ,'com'};
settings.clkCenter = {settings.gpsClkCenter ,settings.gloClkCenter ,settings.galClkCenter ,settings.bdsClkCenter ,'com'};

settings.orbitInterpMethod = 'lagrange';

%% Constants
settings.c           = 299792458;  % Speed of light in m/s
settings.REarth      = 6371.008e3; % Mean Earth radius in m
%% Primary directory
baseDir = preciseProdDir;
tempDir = [baseDir '\Temp Dir\'];

settings.miceDir = miceDir;
settings.tempDir = tempDir;
%% Raw data directories
% GPS----------------------------------------------------------------------
settings.preciseProdDir = ([baseDir 'Raw Data\Precise Products\']);

settings.waasNstbDir = ('D:\WAAS Data\');
settings.outputs = 'D:\GNSS Data\Outputs\Eph Comp and Hists\';
settings.outputs = ([baseDir 'Outputs\Eph Comp and Hists\']);
% RINEX GPS CNAV message files
settings.rnxCnavDir  = ([baseDir 'Raw Data\CNAV_BRDC\']);

% NGA APC Precise Ephemeris
settings.ngaApcDir   = ([baseDir 'Raw Data\Precise_eph\NGA_APC\']);
settings.igsFinDir   = ([baseDir 'Raw Data\Precise_eph\IGS_FINAL\']);

% High-rate precise clock
settings.gpsClkDir   = ([baseDir 'Raw Data\Precise_eph\']);
settings.gloClkDir   = ([baseDir 'Raw Data\Precise_eph\']);
settings.mgxClkDir   = ([baseDir 'Raw Data\Precise_eph\']);

% GLONASS------------------------------------------------------------------
% IGS Precise Ephemeris
settings.gloIgsDir  = ([baseDir 'Raw Data\Precise_eph\GLONASS\']);

% Receiver tropo solutions
settings.igsTropo = ([baseDir 'Raw Data\tropo\']);

% IAC Bulletins------------------------------------------------------------
settings.gloBulDir   = ([baseDir 'Raw Data\GLONASS_BULLETINS\']);

% GALILEO/BEIDOU-----------------------------------------------------------
% MGEX combined broadcast ephemerides
settings.rnxMgexNavDir = ([baseDir 'Raw Data\RINEX_NAV_MGEX\']);

% MGEX Precise ephemerides

% IGS Differential Code Bias data
settings.dcbMgexDir = [baseDir 'Raw Data\DCB\'];

% MGEX Differential Code Bias data
settings.dcbMgexFtp.destDir      = settings.dcbMgexDir;

% IGS Earth Rotation Parameter Data
settings.erpMgexDir = [baseDir 'Raw Data\ERP\'];

% IGS station position solutions
% Receiver positions-------------------------------------------------------
settings.igsStatPos = ([baseDir 'Raw Data\IGS_Stations_Final\']);

% MGEX Precise ephemeris
settings.rnxMgexPephDir = ([baseDir 'Raw Data\Precise_eph\MGEX\']);

% Observation data
settings.rinexObsDir = ([obsDir 'RINEX_OBS\']);
settings.rinexObsMatDir = ([obsDir 'RINEX_OBS_MAT\']);

settings.mgxObsDir =([obsDir 'RINEX_OBS_MGEX\']);
settings.mgxObsMatDir =([obsDir 'RINEX_OBS_MGEX_MAT\']);

settings.mgxHRObsDir = ([obsDir 'RINEX_HR_OBS_MGEX\']);
settings.mgxHRObsMatDir = ([baseDir 'GNSS Data\Obs\RINEX_HR_OBS_MGEX_MAT\']);
settings.mgxHRObsMatDir = ([obsDir 'RINEX_HR_OBS_MGEX_MAT\']);

%% Processed product directories
% SUGL---------------------------------------------------------------------
settings.suglGpsDir = ([baseDir 'Outputs\sugl\sugl_GPS\']);
settings.suglGloDir = ([baseDir 'Outputs\sugl\sugl_GLO\']);
settings.suglGalDir = ([baseDir 'Outputs\sugl\sugl_GAL\']);
settings.suglBdsDir = ([baseDir 'Outputs\sugl\sugl_BDS\']);

% Raw ephemeris comparisons------------------------------------------------
% GPS NGA-COD high rate clock interpolated comparisons
settings.ngaCodCmpDir = ([baseDir 'Outputs\Eph Comp\NGA_COD_CMP\']);
% GPS NGA low rate clock comparisons
settings.ngaCmpDir    = ([baseDir 'Outputs\Eph Comp\NGA_CMP\']);
% Galileo IGS low rate comparisons
settings.mgxGalCmpDir = ([baseDir 'Outputs\Eph Comp\MGX_GAL_CMP\']);
% Beidou IGS low rate comparisons
settings.mgxBdsCmpDir = ([baseDir 'Outputs\Eph Comp\MGX_BDS_CMP\']);


% Histograms---------------------------------------------------------------

% GPS statistics from comparisons
settings.histGpsDir = ([baseDir 'Outputs\histo\histo_GPS\']);
% Galileo statistics from comparisons
settings.histGalDir = ([baseDir 'Outputs\histo\histo_GAL\']);
% Beidou statistics from comparisons
settings.histBdsDir = ([baseDir 'Outputs\histo\histo_BDS\']);

% Receiver observation outputs---------------------------------------------
settings.rxObsDir   = ([baseDir 'Outputs\rxobs\']);

% GLONASS Ephemerides from UNAVCO logs
settings.uraGloDir     = ([baseDir 'Outputs\sugl\GLONASS Eph from T02\']);

% Daily reports
settings.ReportGpsDir = ([baseDir 'Outputs\Daily Reports\GPS Reports\']);
settings.ReportGloDir = ([baseDir 'Outputs\Daily Reports\GLO Reports\']);

%
settings.glonassUploadTimesDir = ([baseDir 'Outputs\misc\GLO Upload Times\']);

% Clock estimation outputs
settings.clockEstOut = ([baseDir 'Misc Outputs\Clock and DCB estimation\']);

end






















