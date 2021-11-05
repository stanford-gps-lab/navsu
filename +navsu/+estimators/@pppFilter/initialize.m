function measId = initialize(obj,obs,corrData,varargin)


p = inputParser;

p.addParameter('measExclude',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
measExclude = res.measExclude;

%% Sort out what is available in the measurements

gnssMeas = navsu.ppp.pullMeasFromList(obs,navsu.internal.MeasEnum.GNSS);
posMeas  = navsu.ppp.pullMeasFromList(obs,navsu.internal.MeasEnum.Position);

if isempty(gnssMeas)
    % Currently, we need a GNSS measurement in order to proceed.
    return;
end

constsUnique = unique(gnssMeas.constInds);


%% Setup state, covariance, and state index mapping
PARAMS = obj.PARAMS;
if isempty(obj.INDS_STATE)
    
    rangeStruc = gnssMeas.range;
    
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
                % If using separate DCBs for each GPS satellite, no need
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

    % Initial guesses at position and velocity
    obj.pos = zeros(3,1);
    
    obj.vel = zeros(3,1);
    
    
    obj.R_b_e = eye(3);
    
    obj.imuBiasStates = [0 0 0 0 0 0]';
    obj.imuBiasStates(1:3) = [-0.3 0.1 -0.4];

    
    obj.clockBias  = zeros(length(constsUnique),1);
    obj.clockDrift =  zeros(length(constsUnique),1);
    
    obj.state = zeros(lastInd,1);
    
    obj.cov = zeros(lastInd,lastInd,1);
end


%%
epoch = gnssMeas.epochs;

% [state,dstate,~,~,covState,covdState,satsUsed] = navsu.ppp.lsSolGnss(gnssMeas,corrData,PARAMS,'pos0',obj.pos);

% Least squares solution, if completed, will populate the position,
% velocity, clock bias, and clock drift states as well as their
% covariances using a simple least squares code phase solution. 
[complete,measId] = leastSquaresSol(obj,epoch,obs,corrData,'measExclude',measExclude);


if ~complete 
    % Unable to initialize given the measurements available
    return;
end


% Initialize attitude - euler angles with body frame

% xi  = -obj.vel./norm(obj.vel);
% z0i = obj.pos./norm(obj.pos);
% yi  = -cross(xi,z0i);
% zi = cross(xi,yi);
obj.R_b_e = navsu.ppp.posVel2Rbe(obj.pos,obj.vel);

obj.imuBiasStates = zeros(6,1);

obj.epochLastInertialUpdate = epoch;
obj.epochLastGnssUpdate     = epoch;

% Check for what satellites have been used 
measIdGnss = measId([measId.TypeID ]== navsu.internal.MeasEnum.GNSS);
satsUsed = unique([cat(1,measIdGnss.prn) cat(1,measIdGnss.const)],'rows');

obj.allSatsSeen = satsUsed;

%% Initialize covariance and state
% Populate the remaining fields in the covariance
cov = obj.cov;

% Attitude
cov(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE)       = eye(3) * PARAMS.SIGMA0.ATTITUDE^2;

% IMU biases
cov(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.ACC_BIAS)       = eye(3) * PARAMS.SIGMA0.ACC_BIAS^2;
cov(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.W_BIAS)           = eye(3) * PARAMS.SIGMA0.W_BIAS^2;

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

%% Initialize phase windup if we have PRN and constInds

prnAll = gnssMeas.PRN(:);
constIndAll = gnssMeas.constInds(:);

if ~isempty(prnAll) && isempty(obj.phWind.phaseOffset)
    nSv = size(prnAll,1);
    obj.phWind.phaseOffset = zeros(size(prnAll,1), 1);
    obj.phWind.PrnConstInd = [prnAll constIndAll];
end

%% Initialization completed successfully!
complete = true;

obj.initialized = true;
end















