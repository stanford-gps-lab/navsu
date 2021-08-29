function [complete,measId,extraInputs] = leastSquaresSol(obj,epoch,obs,corrData,varargin)

p = inputParser;

p.addParameter('measExclude',[]);
p.addParameter('extraInputs',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
measExclude = res.measExclude;
extraInputs0 = res.extraInputs;

%%
state = obj.state;
cov = obj.cov;
nState = length(state);

% Initialize
obj.resids      = [];
obj.measRemoved = [];

gnssMeasMaskCode = obj.PARAMS.measMask;
gnssMeasMaskCode({'pr'},:) = repelem({1},1,size(gnssMeasMaskCode,2));
gnssMeasMaskCode({'cp'},:) = repelem({0},1,size(gnssMeasMaskCode,2));
gnssMeasMaskCode({'dopp'},:) = repelem({1},1,size(gnssMeasMaskCode,2));

% Save initial values in case this doesn't work
pos0        = obj.pos;
vel0        = obj.vel;
clockBias0  = obj.clockBias;
clockDrift0 = obj.clockDrift;

%% Loop through each measurement type and add it to the list!
convMetric = Inf;
convThresh = 1e-3;
while convMetric > convThresh
    
    predMeas = [];
    H = zeros(0,nState);
    R = [];
    measId  = [];       % MeasID
    meas    = [];       % actual measurements
    
    for idx = 1:length(obs)
        obsi = obs{idx};
        
        switch obsi.type
            case navsu.internal.MeasEnum.GNSS
                % this should be code only
                obsi = navsu.ppp.measMask(obsi,gnssMeasMaskCode);
                [predMeasi,Hi,Ri,el,az,prnConstInds,measIdi,measi,~,extraInputs] = ...
                    handleGnssMeas(obj,epoch,obsi,corrData,'SimpleModel',true,...
                    'extraInputs',extraInputs0);
                
            case navsu.internal.MeasEnum.Position
                [predMeasi,Hi,Ri,measIdi,measi] = handlePositionMeas(obj,obsi);
                
            case navsu.internal.MeasEnum.Velocity
                [predMeasi,Hi,Ri,measIdi,measi] = handleVelocityMeas(obj,obsi);
            otherwise
                continue;
        end
        [predMeas,H,R,measId,meas] = catMeas(predMeas,predMeasi,H,Hi,R,Ri,measId,measIdi,meas,measi);
    end
    
    %% State vector just includes position, velocity, clock bias, and clock rate
    indsStateEst = [obj.INDS_STATE.POS obj.INDS_STATE.VEL obj.INDS_STATE.CLOCK_BIAS obj.INDS_STATE.CLOCK_DRIFT];
    
    A = H(:,indsStateEst);
    
    % Check if the system is observable
    if rank(A) < size(A,2)
        break;
    end
    
    resids = meas-predMeas;
    
    if any(isnan(R(:)))
        W = diag(ones(size(A,1),1));
    else
        W = inv(R);
    end
    
    % If there are measurements to exclude, do so now.
    if ~isempty(measExclude)

        measMask = false(size(A,1),1);
        for mdx = 1:length(measExclude)
            measMask = measMask | matches(measExclude(mdx),measId);
        end
        
        % Pull these values out from stuff :)
        A(measMask,:) = [];
        resids(measMask) = [];
        W(measMask,:) = [];
        W(:,measMask) = [];
    end
    
    x =(A'*W*A)\A'*W* resids;
    
    % Update the states -
    obj.pos = obj.pos-x(1:3);
    obj.vel = obj.vel-x(4:6);
    obj.clockBias = obj.clockBias+x(7:(7+length(obj.INDS_STATE.CLOCK_BIAS)-1));
    obj.clockDrift = obj.clockDrift+x((end-length(obj.INDS_STATE.CLOCK_DRIFT)+1):end);
    
    % Convergence metric is just over position and velocity states
    convMetric = norm(x(1:6));
    
end

if convMetric < convThresh
    % the estimation process was a success- keep everything
    
    %% Build the covariance
    covState = inv(A'*W*A);
    
    covFull = obj.cov;
    covFull(indsStateEst,obj.INDS_STATE.POS) = covState(:,1:3);
    covFull(indsStateEst,obj.INDS_STATE.VEL) = covState(:,4:6);
    covFull(indsStateEst,obj.INDS_STATE.CLOCK_BIAS) = covState(:,7:(7+length(obj.INDS_STATE.CLOCK_BIAS)-1));
    covFull(indsStateEst,obj.INDS_STATE.CLOCK_DRIFT) = covState(:,(end-length(obj.INDS_STATE.CLOCK_DRIFT)+1):end);
    
    obj.cov = covFull;
    
    complete = true;
else
    % Return values to how they were before the estimation started
    obj.pos        = pos0;
    obj.vel        = vel0;
    obj.clockBias  = clockBias0;
    obj.clockDrift = clockDrift0;
    
    
    complete = false;
end

end



function [predMeas,H,R,measId,meas] = catMeas(predMeas,predMeasi,H,Hi,R,Ri,measId,measIdi,meas,measi);


% concatenate the measurement information!
% Add everything
predMeas = [predMeas; predMeasi];
H         = [H; Hi];
measId    = [measId; measIdi];
meas      = [meas; measi];

R2 = zeros(size(R,1)+size(Ri,1),size(R,1)+size(Ri,1));
R2(1:size(R,1),1:size(R,1)) = R;
R2((end-(size(Ri,1)-1)):end,(end-(size(Ri,1)-1)):end) = Ri;
R = R2;


end





