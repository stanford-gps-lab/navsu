function [measId, extraInputs] = update(obj,epoch,obs,corrData,varargin)


p = inputParser;

p.addParameter('measExclude',[]);
p.addParameter('extraInputs',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
measExclude = res.measExclude;
extraInputs = res.extraInputs;

%%
gnssMeas = navsu.ppp.pullMeasFromList(obs, navsu.internal.MeasEnum.GNSS);

if isempty(gnssMeas)
    % Currently, we need a GNSS measurement in order to proceed.
    return;
end


%%
if isempty(obj.INDS_STATE)
    constsUnique = unique(gnssMeas.constInds);
    
    PARAMS = obj.PARAMS;
    
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
[complete,measId,extraInputs] = obj.leastSquaresSol(epoch,obs,corrData,...
    'measExclude',measExclude,'extraInputs',extraInputs);

obj.initialized = 2;

obj.measRemoved.id = [];
obj.measRemoved.reason = [];
obj.measRemoved.epoch  = epoch*ones(size([]));

end















