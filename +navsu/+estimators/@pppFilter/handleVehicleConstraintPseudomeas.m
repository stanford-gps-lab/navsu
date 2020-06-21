function [predMeas,H,R,measId,meas] = handleVehicleConstraintPseudomeas(obj)
nState = size(obj.state,1);

velRx = obj.vel;
R_b_e = obj.R_b_e;

%% Initialize outputs
predMeas = zeros(0,1);


%% Previous values for integration :)
if isempty(obj.wheelInfo.R_b_e)
    % the necessary values have not been saved yet- maybe try again next
    % time
    
    H = zeros(0,nState);
    R = zeros(0,1);
    meas = zeros(0,1);
    measId = [];
    
    obj.wheelInfo.vup_int = 0;
    obj.wheelInfo.Hv1_int = zeros(1,3);
    obj.wheelInfo.Hv2_int = zeros(1,3);
    obj.wheelInfo.vcr_int = 0;
    obj.wheelInfo.Hc1_int = zeros(1,3);
    obj.wheelInfo.Hc2_int = zeros(1,3);
    obj.wheelInfo.dt_int   = 0;
    
    return;
end

%%
predMeas = zeros(2,1);
H = zeros(2,nState);
R = zeros(2,2);
meas = zeros(2,1);

%% Vertical velocity constraint
% Predicted "measurement" is the integrated vertical velocity over the time
% increment
vel_up_pred = obj.wheelInfo.vup_int;
meas(1) = -vel_up_pred;
predMeas(1,1) = 0;
measId(1,1) = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.NoSlipVertical);

% Hi = zeros(1,nState);
H(1,obj.INDS_STATE.ATTITUDE) = -obj.wheelInfo.Hv1_int;
H(1,obj.INDS_STATE.VEL)      = -obj.wheelInfo.Hv2_int;

R(1,1) = 0.05^2;

%% cross track velocity constraint
% Predicted "measurement" is the integrated vertical velocity over the time
% increment
vel_cr_pred = obj.wheelInfo.vcr_int;
meas(2) = -vel_cr_pred;
predMeas(2,1) = 0;
measId(2,1) = navsu.internal.MeasIdWheels(1,navsu.internal.MeasEnum.NoSlipCross);

% Hi = zeros(1,nState);
H(2,obj.INDS_STATE.ATTITUDE) = -obj.wheelInfo.Hc1_int;
H(2,obj.INDS_STATE.VEL)      = -obj.wheelInfo.Hc2_int;

R(2,2) = 0.05^2;


%% Zero out the integrals
obj.wheelInfo.vup_int = 0;
obj.wheelInfo.Hv1_int = zeros(1,3);
obj.wheelInfo.Hv2_int = zeros(1,3);
obj.wheelInfo.vcr_int = 0;
obj.wheelInfo.Hc1_int = zeros(1,3);
obj.wheelInfo.Hc2_int = zeros(1,3);
obj.wheelInfo.dt_int   = 0;

end