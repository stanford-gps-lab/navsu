function [predMeas,H,R,measIdi,measi] = handleVehicleConstraintPseudomeas(obj)

nState = size(obj.state,1);

velRx = obj.vel;
R_b_e = obj.R_b_e;

%% Initialize outputs
predMeas = zeros(2,1);
H = zeros(2,nState);
R = zeros(2,2);
measi = zeros(2,1);

% Lever arm from IMU to rear left wheel
lArm = [0.0250, -0.0690, -0.1550]'+[0 1.6256/2 -.2286]';

w = obj.lastGyroMeas;

velWheel = R_b_e'*velRx-navsu.geo.crossProdMatrix(w)*lArm;

%% No vertical velocity constraint (vehicle only moves forward)
pseudoMeasi = 0;
predMeas(1) = [0 0 1]*velWheel;

R(1,1) = .5^2;

H(1,obj.INDS_STATE.VEL) = -[0 0 1]*R_b_e';

measi(1) = pseudoMeasi;

%% No slipping to the right or left constraint
pseudoMeasi = 0;
predMeas(2) = [0 1 0]*velWheel;

R(2,2) = .5^2;

H(2,obj.INDS_STATE.VEL) = -[0 1 0]*R_b_e';

measi(2) = pseudoMeasi;

measIdi = navsu.internal.MeasIdVehicleConstraint([1 1]',[navsu.internal.MeasEnum.NoSlipVertical; navsu.internal.MeasEnum.NoSlipCross]);

end