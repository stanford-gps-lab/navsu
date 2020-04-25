function update(obj,epoch,obs,corrData)

PARAMS = obj.PARAMS;

%% Did the filter initialize successfully
% complete = false;
% just position, velocity, clock bias, clock drift
pos = nan(3,1);
vel = nan(3,1);
clockBias = [];
clockDrift = [];

%%
gnssMeas = navsu.ppp.pullMeasFromList(obs,navsu.internal.MeasEnum.GNSS);


%% If any values weren't provided directly, need to just try some least squares
% Did the filter initialize successfully
[state,dstate,~,~,covState,covdState,satsUsed] = navsu.ppp.lsSolGnss(gnssMeas,corrData,PARAMS);

pos = state(1:3);
vel = dstate(1:3);
clockBias = state(4:end);
clockDrift = dstate(4:end);

constsUnique = unique(gnssMeas.constInds);

%% Populate various fields based on what we have
obj.pos = pos;

obj.vel = vel;

obj.clockBias  = clockBias;
obj.clockDrift = clockDrift;

obj.allSatsSeen = sortrows(satsUsed,2);

INDS_STATE = [];
INDS_STATE.ATTITUDE = 1:3;
INDS_STATE.VEL = 4:6;
INDS_STATE.POS = 7:9;
INDS_STATE.ACC_BIAS = 10:12;
INDS_STATE.W_BIAS = 13:15;

% Add the correct number of clock biases and drifts
lastInd = 15;
INDS_STATE.CLOCK_BIAS = (lastInd+1):(lastInd+length(constsUnique));
INDS_STATE.CLOCK_BIAS_CONSTS = constsUnique;
lastInd = lastInd+length(constsUnique);
INDS_STATE.CLOCK_DRIFT = (lastInd+1):(lastInd+length(constsUnique));
INDS_STATE.CLOCK_DRIFT_CONSTS = constsUnique;
lastInd = lastInd+length(constsUnique);

%% Initialize covariance and state
% Attitude
cov(INDS_STATE.ATTITUDE,INDS_STATE.ATTITUDE)       = eye(3) * PARAMS.SIGMA0.ATTITUDE^2;

% Velocity
cov(INDS_STATE.VEL,INDS_STATE.VEL)             = covdState(1:3,1:3);

% Position
cov(INDS_STATE.POS,INDS_STATE.POS)             = covState(1:3,1:3);

% Clock bias
cov(INDS_STATE.CLOCK_BIAS,INDS_STATE.CLOCK_BIAS)   = covState(4:end,4:end);

% Clock rate
cov(INDS_STATE.CLOCK_DRIFT,INDS_STATE.CLOCK_DRIFT) = covdState(4:end,4:end);

obj.cov = cov;

end















