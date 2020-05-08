function [predMeas,H,R,measIdi,measi] = handleVehicleConstraintPseudomeas(obj)

nState = size(obj.state,1);

velRx = obj.vel;
R_b_e = obj.R_b_e;

%% Initialize outputs
predMeas = zeros(2,1);
H = zeros(2,nState);
R = zeros(2,2);
measi = zeros(2,1);

%% No vertical velocity constraint (vehicle only moves forward)
pseudoMeasi = 0;
predMeas(1) = [0 0 1]*R_b_e'*velRx;

R(1,1) = 1^2;

H(1,obj.INDS_STATE.VEL) = -[0 0 1]*R_b_e';

measi(1) = pseudoMeasi;

%% No slipping to the right or left constraint
pseudoMeasi = 0;
predMeas(2) = [0 1 0]*R_b_e'*velRx;

R(2,2) = 1^2;

H(2,obj.INDS_STATE.VEL) = -[0 1 0]*R_b_e';

measi(2) = pseudoMeasi;

measIdi = navsu.internal.MeasIdVehicleConstraint([1 1]',[navsu.internal.MeasEnum.NoSlipVertical; navsu.internal.MeasEnum.NoSlipCross]);

end