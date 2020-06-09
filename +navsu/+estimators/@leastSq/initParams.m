function PARAMS = initParams(obj)

PARAMS = [];

% Flag for stationary scenarios
PARAMS.stationaryMode = false;

% Sub-structure relating to what measurements are to be used
PARAMS.measUse = struct(...
    'dfOnly',            false,...         % true = dual frequency meas only
    'SIG1_THRESH',         20,...            % SNR threshold for L1
    'SIG2_THRESH',         20,...            % SNR threshold for L2
    'SIG3_THRESH',         20,...            % SNR threshold for L2
    'gFreeSlipThresh',   0.05,...          % Threshold for cycle slip- geometry free [m]
    'slipDetector',      'GFREE',...   % Cycle slip detector ('RX_OUTPUT','GFREE')
    'noVertVel',         false);          % 0 vertical velocity constraint


PARAMS.measUse.excludeThresh.(char(navsu.internal.MeasEnum.GNSS)) = 5;
PARAMS.measUse.excludeThreshLarge.(char(navsu.internal.MeasEnum.GNSS)) = 20;

PARAMS.measUse.excludeThresh.(char(navsu.internal.MeasEnum.Position)) = Inf;
PARAMS.measUse.excludeThreshLarge.(char(navsu.internal.MeasEnum.Position)) = Inf;

PARAMS.measUse.excludeThresh.(char(navsu.internal.MeasEnum.Position)) = Inf;
PARAMS.measUse.excludeThreshLarge.(char(navsu.internal.MeasEnum.Position)) = Inf;

PARAMS.measUse.excludeThresh.(char(navsu.internal.MeasEnum.NoSlipCross)) = Inf;
PARAMS.measUse.excludeThreshLarge.(char(navsu.internal.MeasEnum.NoSlipCross)) = Inf;

PARAMS.measUse.excludeThresh.(char(navsu.internal.MeasEnum.NoSlipVertical)) = Inf;
PARAMS.measUse.excludeThreshLarge.(char(navsu.internal.MeasEnum.NoSlipVertical)) = Inf;

% Measurement masking- which measurments to actually use in the
% filter
PARAMS.measMask = table([1 1 1]',[1 0 0]',[1 0 0]',[1 1 0]',[1 1 0]',...
    'VariableNames',{'f1','f2','f3','f1f2','f1f3'},...
    'RowNames',{'pr','cp','dopp'});

% measurement uncertainty
PARAMS.sigMeas = struct(...
    'pr',                1,...             % pseudorange   [m]
    'cp',             0.03,...             % carrier phase [m]
    'dopp',           0.05,...             % doppler       [m/s]
    'iono',           1000,...             % iono sigma- no correction [m]
    'ionoRate',        100);               % iono rate sigma- no correction [m/s]

% What states to include in the filter and probably some settings
% here
PARAMS.states = struct(...
    'attitudeMode','EULER',...        % Attitude representation
    'RX_DCB',       true,...          % Receiver DCB for multi-frequency
    'trop',         true,...          % Tropospheric state
    'iono',         true,...          % Iono state for single frequency meas
    'ionoMode',     'L1DELAYSTATE',...% L1DELAYSTATE estimates delay at L1 per LOS
    'RX_DCB_GLO',   true,...          % Separate DCB state for each GLONASS code measurement
    'RX_DCB_GPS',   false,...         % Separate DCB state for each GLONASS code measurement
    'MP_CODE',      false,...          % Code phase multipath
    'MP_CARR',      false);            % Carrier phase multipath

% Tropospheric model
PARAMS.tropModel = 'UNB3';

% Satellite DCB use flag
PARAMS.dcbUse = true;

% Body frame IMU lever arm - IMU to GNSS antenna
PARAMS.IMU_ARM = [0 0 0]';

% Elevation mask angle
PARAMS.elMask = 7.5; % Degrees

% Speed of light in meters
PARAMS.c = 299792458;

% Initial uncertainties
PARAMS.SIGMA0 = struct(...
    'POS',      5,...             % Position
    'VEL',      0.1,...           % Velocity
    'ACC',      [2.5 2.5 0.25],...% this is actually not used anymore
    'ATTITUDE', 3*pi/180,...     % Attitude
    'ACC_BIAS', 1e-2*70e-3,...    % Accelerometer bias
    'W_BIAS',   1e-2*1*pi/180,... % Gyro bias
    'ATT_RATE', 1e-0*2*pi/180,... % Attitude rate
    'ACC_SCALE',0.04,...          % Accelerometer scale factor
    'W_SCALE',  0.01,...          % Gyro scale
    'TROP',     0.05,...          % Tropospheric delta state
    'AMB',      100,...             % Carrier phase ambiguity
    'RX_DCB',   5,...             % Receiver DCB
    'L1_IONO',  1,...             % Slant L1 iono delay
    'RX_DCB_GLO', 1,...           % Separate DCB state for each GLONASS code measurement
    'RX_DCB_GPS', 1,...           % Separate DCB state for each GPS code measurement
    'MP_CODE',    2,...           % Code phase multipath
    'MP_CARR',    0.03);          % Carrier phase multipath

% Process noise (need to be squared later :) )
PARAMS.Q = struct(...
    'POS',             1,...                 % Position
    'VEL',             0.5,...               % Velocity
    'ACC',             3.2,...               % Acceleration
    'ATTTITUDE',       [20 20 20]*pi/180,... % Attitude
    'ATT_RATE',        [40 40 40]*pi/180,... % Atttitude rate
    'W_BIAS',          2e-5*1,...%/100,...   % Gyro bias
    'ACC_BIAS',        1e-3*1,...%/100,...             % Accelerometer bias
    'ACC_SCALE',       1e-5,...              % Accelerometer scale
    'W_SCALE',         1e-5,...              % Gyro scale
    'RXB',             10,...                % Receiver clock bias
    'gyro_noise_PSD',  0.0015,...%*100,...   % Gyro noise
    'accel_noise_PSD', 0.005,...%*1000,...   % Accelerometer noise
    'RX_DCB',          0,...                 % Receiver DCB
    'TROP',            0.002/60,...          % Tropospheric delta state
    'AMB',             0,...                 % Carrier phase ambiguity
    'L1_IONO',         0.03,...              % Slant L1 iono delay
    'RX_DCB_GLO',      0,...                 % Separate DCB state for each GLONASS code measurement
    'RX_DCB_GPS',      0,...                 % Separate DCB state for each GLONASS code measurement  \
    'MP_CODE',         0.25,....             % Code phase multipath
    'MP_CARR',         0.0000);              % Carrier phase multipath

PARAMS.other = struct(...
    'TAU_MP_CODE',    100);                  % Time constant for code multipath

% Related to solution separation
PARAMS.solSep = struct(...
    'subsetGroupSize', 1,...
    'Psat', 1e-5,...
    'nOutSubset',1,...
    'faultResetSubsets',false,...
    'nMaxSubset',100,...
    'pl_ss_exact_optfa',0,...
    'pl_ss_approx1_optfa',1,...
    'pl_ss_approx2_optfa',0,...
    'pl_ss_exact',0,...
    'pl_ss_approx1',0,...
    'pl_ss_approx2',0,...
    'fdeMode','reinit');  % 'reinit' (just start over when a fault is detected) or 'background'

end