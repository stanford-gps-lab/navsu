function complete = initialize(obj,PARAMS,corrData,varargin)

p = inputParser;

p.addParameter('pos',[]);
p.addParameter('vel',[]);
p.addParameter('R_b_e',[]);
p.addParameter('imuBiasStates',[]);
p.addParameter('constsUnique',1);
p.addParameter('clockBias',[]);
p.addParameter('clockDrift',[]);
p.addParameter('PRN',[]);
p.addParameter('constInds',[]);
p.addParameter('range',[]);
p.addParameter('epoch',[]);
p.addParameter('gnssMeas',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
pos           = res.pos;
vel           = res.vel;
R_b_e         = res.R_b_e;
imuBiasStates = res.imuBiasStates;
constsUnique  = res.constsUnique;
clockDrift    = res.clockDrift;
clockBias     = res.clockBias;
PRN           = res.PRN(:);
constInds     = res.constInds(:);
rangeStruc    = res.range;
gnssMeas      = res.gnssMeas;
epoch         = res.epoch;

rangeStruc = gnssMeas.range;
%% Did the filter initialize successfully
% complete = false;

%% If any values weren't provided directly, need to just try some least squares
% Did the filter initialize successfully
complete = ~any(cellfun(@isempty,{pos vel R_b_e clockBias clockDrift}));

covState = [];
covdState = [];
satsUsed = [];
if ~complete && ~isempty(gnssMeas)
    [state,dstate,~,~,covState,covdState,satsUsed] = navsu.ppp.lsSolGnss(gnssMeas,corrData,PARAMS);
    
    if ~any(state)
        % didn't get any outputs here- just escape
        return;
    end
    
    pos0 = state(1:3);
    vel0 = dstate(1:3);
    epoch0 = gnssMeas.epochs;
    clockBias0 = state(4:end);
    clockDrift0 = dstate(4:end);
    
    % Initialize attitude - euler angles with body frame
    xi  = -vel0./norm(vel0);
    z0i = pos0./norm(pos0);
    yi  = -cross(xi,z0i);
    zi = cross(xi,yi);
    R_b_e0 = [xi yi zi];
    
    if any(isnan(R_b_e0(:)))
        R_b_e0 = eye(3);
    end
    
    pos0 = pos0-R_b_e0*PARAMS.IMU_ARM;
    
    % Initialize bias estimates
    imuBiasStates = zeros(6,1); % accelerometer(1:3) and gyro(1:3)
    imuBiasStates(1:3) = [-0.3 0.1 -0.4];
    
    % replace all of the values that didn't have one previously
    if isempty(pos)
        pos = pos0;
    end
    if isempty(vel)
        vel = vel0;
    end
    if isempty(R_b_e)
        R_b_e = R_b_e0;
    end
    if isempty(clockBias)
        clockBias = clockBias0;
    end
    if isempty(clockDrift)
        clockDrift = clockDrift0;
    end
    epoch = epoch0;
    
    constsUnique = unique(gnssMeas.constInds);
    
    PRN = gnssMeas.PRN(:);
    constInds = gnssMeas.constInds(:);
end

%% Populate various fields based on what we have
obj.pos = pos;

% if isempty(obj.vel)
obj.vel = vel;
% end


obj.R_b_e = R_b_e;

obj.imuBiasStates = imuBiasStates;

obj.clockBias  = clockBias;
obj.clockDrift = clockDrift;

obj.epochLastInertialUpdate = epoch;
obj.epochLastGnssUpdate     = epoch;

obj.allSatsSeen = sortrows(satsUsed,2);

%% Setup state, covariance, and state index mapping
if isempty(obj.INDS_STATE)
    obj.INDS_STATE = [];
    obj.INDS_STATE.ATTITUDE = 1:3;
    obj.INDS_STATE.VEL = 4:6;
    obj.INDS_STATE.POS = 7:9;
    obj.INDS_STATE.ACC_BIAS = 10:12;
    obj.INDS_STATE.W_BIAS = 13:15;
    
    % Add the correct number of clock biases and drifts
    lastInd = 15;
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
    
    if PARAMS.states.RX_DCB
        % For each constellation and signal (excluding one reference), we need
        % a receiver differential code bias state
        if isempty(rangeStruc)
            error('Ranging measuremnets required for RX DCB state initialization')
        end
        
        sigConstInds = unique([rangeStruc.sig(:) rangeStruc.constInds(:)],'rows');
        constUn = unique(rangeStruc.constInds(:));
        for idx = 1:length(constUn)
            if constUn(idx) == 1 && PARAMS.states.RX_DCB_GPS
                % If using separate DCBs for each GLONASS satellite, no need
                % for a signal specific one
                continue;
            end
            
            if constUn(idx) == 2 && PARAMS.states.RX_DCB_GLO
                % If using separate DCBs for each GLONASS satellite, no need
                % for a signal specific one
                continue;
            end
            % Remove a reference signal from each constellation- just using the largest one for now.
            % this should change based on whereever measurement masking happens
            refSigi = max(sigConstInds(sigConstInds(:,2) == constUn(idx),1));
            
            sigConstInds(ismember(sigConstInds,[refSigi constUn(idx)],'rows'),:) = [];
        end
        
        % Add the state indices
        nStateDcb = size(sigConstInds,1);
        
        obj.INDS_STATE.RX_DCB.INDS = lastInd+[1:nStateDcb];
        obj.INDS_STATE.RX_DCB.sig  = sigConstInds(:,1)';
        obj.INDS_STATE.RX_DCB.constInds = sigConstInds(:,2)';
        
        lastInd = lastInd+nStateDcb;
    else
        obj.INDS_STATE.RX_DCB.INDS = [];
        obj.INDS_STATE.RX_DCB.sig  = [];
        obj.INDS_STATE.RX_DCB.constInds = [];
    end
    
    obj.INDS_STATE.FLEX_STATE_MIN = lastInd+1;
    obj.INDS_STATE.FLEX_STATES = zeros(0,1);
    obj.INDS_STATE.FLEX_STATES_INFO = zeros(0,4); % PRN | CONST | TYPE | SIG INDICATOR
    
end

%% Initialize covariance and state
if isempty(obj.state)
    state = zeros(lastInd,1);
else
    state = zeros(size(obj.state));
    %     state = obj.state;
end

if isempty(obj.cov)
    cov = zeros(lastInd);
else
    cov = obj.cov;
    cov = diag(1000*ones(size(obj.cov,1),1));
end

% Attitude
cov(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE)       = eye(3) * PARAMS.SIGMA0.ATTITUDE^2;

% Velocity
if isempty(covdState)
    cov(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL)             = eye(3) * PARAMS.SIGMA0.VEL^2;
else
    cov(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL)             = covdState(1:3,1:3);
end

% Position
if isempty(covState)
    cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS)             = eye(3) * PARAMS.SIGMA0.POS^2;
else
    cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS)             = covState(1:3,1:3);
end

% IMU biases
cov(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.ACC_BIAS)       = eye(3) * PARAMS.SIGMA0.ACC_BIAS^2;
cov(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.W_BIAS)           = eye(3) * PARAMS.SIGMA0.W_BIAS^2;

% Clock bias
if isempty(covState)
    cov(obj.INDS_STATE.CLOCK_BIAS,obj.INDS_STATE.CLOCK_BIAS)   = eye(length(obj.INDS_STATE.CLOCK_BIAS))*100^2;
else
    cov(obj.INDS_STATE.CLOCK_BIAS,obj.INDS_STATE.CLOCK_BIAS)   = covState(4:end,4:end);
end

% Clock rate
if isempty(covdState)
    cov(obj.INDS_STATE.CLOCK_DRIFT,obj.INDS_STATE.CLOCK_DRIFT) = eye(length(obj.INDS_STATE.CLOCK_DRIFT))*100^2;
else
    cov(obj.INDS_STATE.CLOCK_DRIFT,obj.INDS_STATE.CLOCK_DRIFT) = covdState(4:end,4:end);
end

% Tropo
cov(obj.INDS_STATE.TROP,obj.INDS_STATE.TROP)               = PARAMS.SIGMA0.TROP^2;
% DCBs
cov(obj.INDS_STATE.RX_DCB.INDS,obj.INDS_STATE.RX_DCB.INDS) = eye(length(obj.INDS_STATE.RX_DCB.INDS))*PARAMS.SIGMA0.RX_DCB^2;

% There may actually be additional states if this is being initialized
% from another filter.
if ~isempty(obj.INDS_STATE.FLEX_STATES)
    state2 = nan(size(obj.state ));
    covDiag2 = nan(size(obj.state));
    stateTypes = {'cp'; ... % 1
        'L1DELAYSTATE'; ... % 2
        'MP_CODE';...       % 3
        'MP_CARR';...       % 4
        'RX_DCB_GLO'};      % 5
    
    % Go through and initialize each state haha
    for idx = 1:size(obj.INDS_STATE.FLEX_STATES)
        stateTypei = stateTypes{obj.INDS_STATE.FLEX_STATES_INFO(idx,3)};
        
        [stateAdd,covAdd] = initStateCov(stateTypei,obj.INDS_STATE.FLEX_STATES_INFO(idx,:),PARAMS,gnssMeas);
        
        state(obj.INDS_STATE.FLEX_STATES(idx)) = stateAdd;
        cov(obj.INDS_STATE.FLEX_STATES(idx),obj.INDS_STATE.FLEX_STATES(idx)) = covAdd;
        
    end
end

obj.cov = cov;

obj.state = state;

%% Initialize phase windup if we have PRN and constInds
if ~isempty(PRN) && ~isempty(constInds) && isempty(obj.phWind.phaseOffset)
    nSv = length(PRN);
    obj.phWind.phaseOffset = zeros(size(PRN));
    obj.phWind.PrnConstInd = [PRN constInds];
end

%% Initialization completed successfully!
complete = true;


obj.initialized = true;
end















