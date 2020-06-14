function [predMeas,H,R,measId,meas] = handleWheelsMeas(obj,wheelMeas)

nState = size(obj.state,1);

velRx = obj.vel;
R_b_e = obj.R_b_e;

%% Initialize outputs
predMeas = zeros(0,1);

% % Lever arm from IMU to rear left wheel
lArm = [0.0250, -0.0690, -0.1550]'+[0 1.6256/2 -.2286]';
%
% %%
% 30.5 inch wheel
cWheel = 30.5*2.54/100*pi;
scale0 = cWheel/1000;

%% Scale factors
sFactors = obj.state(obj.INDS_STATE.WHEELS);
sfl = sFactors(1);  %front left
sfr = sFactors(2);  %front right
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
    obj.wheelInfo.vup_int = 0;
    obj.wheelInfo.Hv1_int = zeros(1,3);
    obj.wheelInfo.Hv2_int = zeros(1,3);
    obj.wheelInfo.vcr_int = 0;
    obj.wheelInfo.Hc1_int = zeros(1,3);
    obj.wheelInfo.Hc2_int = zeros(1,3);
    
    return;
end

%%
predMeas = zeros(2,1);
H = zeros(2,nState);
R = zeros(2,2);
meas = zeros(2,1);

R_b_e_prev = obj.wheelInfo.R_b_e;
vel_prev   = obj.wheelInfo.vel;
w_prev     = obj.wheelInfo.w;

% current values
R_b_e = obj.R_b_e;
vel   = obj.vel;
w     = obj.lastGyroMeas;

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

R(1,1) = 0.005^2;

%% Handle the rear right wheel measurement
vel_rr_pred   = obj.wheelInfo.vrr_int;
meas(2)       = scale0*wheelMeas.SpeedRearRight*(1-srr)-vel_rr_pred;
predMeas(2,1) = 0;
measId(2,1)   = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.SpeedRearRight);


H(2,obj.INDS_STATE.ATTITUDE)  = -obj.wheelInfo.H11_int;
H(2,obj.INDS_STATE.VEL)       = -obj.wheelInfo.H12_int;
H(2,obj.INDS_STATE.WHEELS(4)) = vel_rr_pred;

R(2,2) = 0.005^2;

%% Vertical velocity constraint
% Predicted "measurement" is the integrated vertical velocity over the time
% increment
vel_up_pred = obj.wheelInfo.vup_int;
meas(3) = -vel_up_pred;
predMeas(3,1) = 0;
measId(3,1) = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.NoSlipVertical);

% Hi = zeros(1,nState);
H(3,obj.INDS_STATE.ATTITUDE) = -obj.wheelInfo.Hv1_int;
H(3,obj.INDS_STATE.VEL)      = -obj.wheelInfo.Hv2_int;

R(3,3) = 0.1^2;

%% cross track velocity constraint
% Predicted "measurement" is the integrated vertical velocity over the time
% increment
vel_cr_pred = obj.wheelInfo.vcr_int;
meas(4) = -vel_cr_pred;
predMeas(4,1) = 0;
measId(4,1) = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.NoSlipCross);

% Hi = zeros(1,nState);
H(4,obj.INDS_STATE.ATTITUDE) = -obj.wheelInfo.Hc1_int;
H(4,obj.INDS_STATE.VEL)      = -obj.wheelInfo.Hc2_int;

R(4,4) = 0.1^2;



%% Zero out the integrals
obj.wheelInfo.vrl_int = 0;
obj.wheelInfo.vrr_int = 0;
obj.wheelInfo.H11_int = zeros(1,3);
obj.wheelInfo.H12_int = zeros(1,3);
obj.wheelInfo.vup_int = 0;
obj.wheelInfo.Hv1_int = zeros(1,3);
obj.wheelInfo.Hv2_int = zeros(1,3);
obj.wheelInfo.vcr_int = 0;
obj.wheelInfo.Hc1_int = zeros(1,3);
obj.wheelInfo.Hc2_int = zeros(1,3);

if 0 
l_rl = [0.0250, -0.0690, -0.1550]'+[0 0.6256/2 -.2286]';

% this is currently just using a super stupid integration using the value
% at the previous time step lol
vel_rl_prev = [-1 0 0]*(R_b_e_prev'*vel_prev-navsu.geo.crossProdMatrix(w_prev)*l_rl);
vel_rl_curr = [-1 0 0]*(R_b_e'*vel-navsu.geo.crossProdMatrix(w)*l_rl);

vel_rl_pred = 1/2*(vel_rl_prev+vel_rl_curr);


meas(1) = scale0*wheelMeas.SpeedRearLeft*(1-srl)-vel_rl_pred;
predMeas(1) = 0;

measId(1) = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.SpeedRearLeft);
% H = zeros(1,nState);

H11 = -1/1*1/2*([-1 0 0]*R_b_e_prev'*navsu.geo.crossProdMatrix(vel_prev)+...
    [-1 0 0]*R_b_e'*navsu.geo.crossProdMatrix(vel));

H12 = -1/1*1/2*([-1 0 0]*R_b_e_prev'+[-1 0 0]*R_b_e');

H(1,obj.INDS_STATE.ATTITUDE) = H11;
H(1,obj.INDS_STATE.VEL) = H12;
H(1,obj.INDS_STATE.WHEELS(3)) = vel_rl_pred;

R(1,1) = 0.5^2;

%% Predict the rear right wheel measurement
l_rl = [0.0250, -0.0690, -0.1550]'+[0 -1.6256/2 -.2286]';

% this is currently just using a super stupid integration using the value
% at the previous time step lol
vel_rr_prev = [-1 0 0]*(R_b_e_prev'*vel_prev-navsu.geo.crossProdMatrix(w_prev)*l_rl);
vel_rr_curr = [-1 0 0]*(R_b_e'*vel-navsu.geo.crossProdMatrix(w)*l_rl);

vel_rr_pred = 1/2*(vel_rr_prev+vel_rr_curr);

meas(2) = scale0*wheelMeas.SpeedRearLeft*(1-srr)-vel_rr_pred;
predMeas(2) = 0;

measId(2,1) = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.SpeedRearRight);
% H = zeros(1,nState);

H11 = -1/1*1/2*([-1 0 0]*R_b_e_prev'*navsu.geo.crossProdMatrix(vel_prev)+...
    [-1 0 0]*R_b_e'*navsu.geo.crossProdMatrix(vel));

H12 = -1/1*1/2*([-1 0 0]*R_b_e_prev'+[-1 0 0]*R_b_e');

H(2,obj.INDS_STATE.ATTITUDE) = H11;
H(2,obj.INDS_STATE.VEL) = H12;
H(2,obj.INDS_STATE.WHEELS(4)) = vel_rr_pred;

R(2,2) = 0.5^2;
end

end