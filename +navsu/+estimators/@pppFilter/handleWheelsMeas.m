function [predMeas,H,R,measId,meas] = handleWheelsMeas(obj,wheelMeas)

nState = size(obj.state,1);

%% Initialize outputs
predMeas = zeros(0,1);

% 30.5 inch wheel
cWheel = (30.5)*2.54/100*pi;
scale0 = cWheel/1000;

%% Scale factors
sFactors = obj.state(obj.INDS_STATE.WHEELS);
srl = sFactors(3);  %rear left
srr = sFactors(4);  %rear right

%% Previous values for integration :)
if isempty(obj.wheelInfo.R_b_e)
    % the necessary values have not been saved yet- maybe try again next
    % time
    
    H = zeros(0,nState);
    R = zeros(0,1);
    meas = zeros(0,1);
    measId = [];
    
    obj.wheelInfo.vrl_int = 0;
    obj.wheelInfo.vrr_int = 0;
    obj.wheelInfo.H11_int = zeros(1,3);
    obj.wheelInfo.H12_int = zeros(1,3);
    obj.wheelInfo.dt_int   = 0;
    
    return;
end

%%
predMeas = zeros(2,1);
H = zeros(2,nState);
R = zeros(2,2);
meas = zeros(2,1);

%% Handle the rear left wheel measurement
% Predicted integrated velocity comes from mechanization output
vel_rl_pred   = obj.wheelInfo.vrl_int;
% Measurement is just innovation
meas(1)       = scale0*wheelMeas.SpeedRearLeft*(1-srl)-vel_rl_pred;
% Predicted measurement is thus 0
predMeas(1,1) = 0;
% Measurement ID
measId(1,1)   = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.SpeedRearLeft);

% Measurement sensitivies were also integrated in mechanization
H(1,obj.INDS_STATE.ATTITUDE)  = -obj.wheelInfo.H11_int;
H(1,obj.INDS_STATE.VEL)       = -obj.wheelInfo.H12_int;
H(1,obj.INDS_STATE.WHEELS(3)) = vel_rl_pred;

R(1,1) = 0.05^2;
R(1,1) = 0.1^2;


%% Handle the rear right wheel measurement
vel_rr_pred   = obj.wheelInfo.vrr_int;
meas(2)       = scale0*wheelMeas.SpeedRearRight*(1-srr)-vel_rr_pred;
predMeas(2,1) = 0;
measId(2,1)   = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.SpeedRearRight);

H(2,obj.INDS_STATE.ATTITUDE)  = -obj.wheelInfo.H11_int;
H(2,obj.INDS_STATE.VEL)       = -obj.wheelInfo.H12_int;
H(2,obj.INDS_STATE.WHEELS(4)) = vel_rr_pred;

R(2,2) = 0.05^2;
R(2,2) = 0.1^2;


%% Zero out the integrals
obj.wheelInfo.vrl_int = 0;
obj.wheelInfo.vrr_int = 0;
obj.wheelInfo.H11_int = zeros(1,3);
obj.wheelInfo.H12_int = zeros(1,3);
obj.wheelInfo.dt_int   = 0;


end