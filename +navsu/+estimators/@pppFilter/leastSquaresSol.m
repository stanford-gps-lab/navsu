function [state,cov] = leastSquaresSol(obj,epoch,obs,corrData)

state = obj.state;
cov = obj.cov;
nState = length(state);

% nState =

% Initialize
% predMeas = [];
% H = zeros(0,nState);
% R = [];
% measId  = [];       % MeasID
% meas    = [];       % actual measurements
% prnConstInds = [];
% el = [];
% az = [];

obj.resids      = [];
obj.measRemoved = [];
gnssMeas = [];

gnssMeasMaskCode = obj.PARAMS.measMask;
gnssMeasMaskCode({'pr'},:) = repelem({1},1,size(gnssMeasMaskCode,2));
gnssMeasMaskCode({'cp'},:) = repelem({0},1,size(gnssMeasMaskCode,2));
gnssMeasMaskCode({'dopp'},:) = repelem({1},1,size(gnssMeasMaskCode,2));

%% Loop through each measurement type and add it to the list!
for jdx = 1:10
    
    predMeas = [];
    H = zeros(0,nState);
    R = [];
    measId  = [];       % MeasID
    meas    = [];       % actual measurements
    prnConstInds = [];
    el = [];
    az = [];
    for idx = 1:length(obs)
        obsi = obs{idx};
        
        switch obsi.type
            case navsu.internal.MeasEnum.GNSS
                % this should be code only
                
                obsi = navsu.ppp.measMask(obsi,gnssMeasMaskCode);
                
                [predMeasi,Hi,Ri,el,az,prnConstInds,measIdi,measi] = ...
                    handleGnssMeas(obj,epoch,obsi,corrData,'SimpleModel',true);
                
                gnssMeas = obsi;
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
    
    % 
    indsStateEst = [obj.INDS_STATE.POS obj.INDS_STATE.VEL obj.INDS_STATE.CLOCK_BIAS obj.INDS_STATE.CLOCK_DRIFT];
    
    A = H(:,indsStateEst);
    
    resids = meas-predMeas;
    
    W = diag(ones(size(A,1),1));
    W=  inv(R);
    
    % Need to check for observability of each state
% %     if rank(A) < size(A,2)
% %         
% %         
% %     end
%     W 
%     resids = 
    
    x =(A'*W*A)\A'*W* resids;
    
    x = x;
    % Update the states lol
    obj.pos = obj.pos+x(1:3);
%     obj.vel = obj.vel+x(4:6);
%     obj.clockBias = obj.clockBias+x(7:(7+length(obj.INDS_STATE.CLOCK_BIAS)-1));
%     obj.clockDrift = obj.clockDrift+x((end-length(obj.INDS_STATE.CLOCK_DRIFT)+1):end);
    
    %     state = state+x;
    covState = inv(A'*W*A);
    
    norm(x)
    
    
    
    
    
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