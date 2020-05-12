function update(obj,epoch,obs,corrData)

PARAMS = obj.PARAMS;

%%
gnssMeas = navsu.ppp.pullMeasFromList(obs,navsu.internal.MeasEnum.GNSS);

if isempty(gnssMeas)
    % Currently, we need a GNSS measurement in order to proceed.
    return;
end

constsUnique = unique(gnssMeas.constInds);

%%
if isempty(obj.INDS_STATE)
    constsUnique = unique(gnssMeas.constInds);
    
    PARAMS = obj.PARAMS;
    rangeStruc = gnssMeas.range;
    
    obj.INDS_STATE = [];
    obj.INDS_STATE.ATTITUDE = 1:3;
    obj.INDS_STATE.VEL = 4:6;
    obj.INDS_STATE.POS = 7:9;
    
    % Add the correct number of clock biases and drifts
    lastInd = 9;
    obj.INDS_STATE.CLOCK_BIAS = (lastInd+1):(lastInd+length(constsUnique));
    obj.INDS_STATE.CLOCK_BIAS_CONSTS = constsUnique;
    lastInd = lastInd+length(constsUnique);
    obj.INDS_STATE.CLOCK_DRIFT = (lastInd+1):(lastInd+length(constsUnique));
    obj.INDS_STATE.CLOCK_DRIFT_CONSTS = constsUnique;
    lastInd = lastInd+length(constsUnique);
    
    if PARAMS.states.trop
        obj.INDS_STATE.TROP = lastInd+1;
        
        lastInd = lastInd+1;
    else
        obj.INDS_STATE.TROP = [];
    end
    
    obj.INDS_STATE.FLEX_STATE_MIN = lastInd+1;
    obj.INDS_STATE.FLEX_STATES = zeros(0,1);
    obj.INDS_STATE.FLEX_STATES_INFO = zeros(0,4); % PRN | CONST | TYPE | SIG INDICATOR
    
     % Initial guesses at position and velocity
    obj.pos = zeros(3,1);
    
    obj.vel = zeros(3,1);
    
    
    obj.R_b_e = eye(3);
    
    obj.imuBiasStates = [0 0 0 0 0 0]';
    
    obj.clockBias  = zeros(length(constsUnique),1);
    obj.clockDrift =  zeros(length(constsUnique),1);
    
    obj.state = zeros(lastInd,1);
    
    obj.cov = zeros(lastInd,lastInd,1);
end

%% Compute the least squares solution
[complete,measIds] = leastSquaresSol(obj,epoch,obs,corrData);


%% If any values weren't provided directly, need to just try some least squares
% Did the filter initialize successfully
% [state,dstate,~,~,covState,covdState,satsUsed] = navsu.ppp.lsSolGnss(gnssMeas,corrData,PARAMS);
% 
% pos = state(1:3);
% vel = dstate(1:3);
% clockBias = state(4:end);
% clockDrift = dstate(4:end);
% 
% 
% %% Populate various fields based on what we have
% obj.pos = pos;
% 
% obj.vel = vel;
% 
% obj.clockBias  = clockBias;
% obj.clockDrift = clockDrift;
% 
% obj.allSatsSeen = sortrows(satsUsed,2);
% 
% % if isempty(obj.INDS_STATE)
% INDS_STATE = [];
% INDS_STATE.ATTITUDE = 1:3;
% INDS_STATE.VEL = 4:6;
% INDS_STATE.POS = 7:9;
% INDS_STATE.ACC_BIAS = 10:12;
% INDS_STATE.W_BIAS = 13:15;
% 
% % Add the correct number of clock biases and drifts
% lastInd = 15;
% INDS_STATE.CLOCK_BIAS = (lastInd+1):(lastInd+length(constsUnique));
% INDS_STATE.CLOCK_BIAS_CONSTS = constsUnique;
% lastInd = lastInd+length(constsUnique);
% INDS_STATE.CLOCK_DRIFT = (lastInd+1):(lastInd+length(constsUnique));
% INDS_STATE.CLOCK_DRIFT_CONSTS = constsUnique;
% lastInd = lastInd+length(constsUnique);
% 
% obj.INDS_STATE = INDS_STATE;
% 
% %% Initialize covariance and state
% % Attitude
% cov(INDS_STATE.ATTITUDE,INDS_STATE.ATTITUDE)       = eye(3) * PARAMS.SIGMA0.ATTITUDE^2;
% 
% % Velocity
% if ~isempty(covdState)
%     cov(INDS_STATE.VEL,INDS_STATE.VEL)             = covdState(1:3,1:3);
% end
% 
% % Position
% cov(INDS_STATE.POS,INDS_STATE.POS)             = covState(1:3,1:3);
% 
% % Clock bias
% cov(INDS_STATE.CLOCK_BIAS,INDS_STATE.CLOCK_BIAS)   = covState(4:end,4:end);
% 
% % Clock rate
% if ~isempty(covdState)
%     cov(INDS_STATE.CLOCK_DRIFT,INDS_STATE.CLOCK_DRIFT) = covdState(4:end,4:end);
% end
% 
% obj.cov = cov;

end















