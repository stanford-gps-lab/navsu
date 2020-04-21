function [predMeas,H,R,measMat] = handleVehicleConstraintPseudomeas(obj)

nState = size(obj.state,1);

velRx = obj.vel;
R_b_e = obj.R_b_e;

%% Initialize outputs
predMeas = zeros(2,1);
H = zeros(2,nState);
R = zeros(2,2);
measMat = zeros(2,6);

%% No vertical velocity constraint (vehicle only moves forward)
pseudoMeasi = 0;
predMeas(1) = [0 0 1]*R_b_e'*velRx;

R(1,1) = 1^2;

H(1,obj.INDS_STATE.VEL) = -[0 0 1]*R_b_e';

measMat(1,:) = [0 0 0 0 pseudoMeasi 4];


%% No slipping to the right or left constraint
pseudoMeasi = 0;
predMeas(2) = [0 1 0]*R_b_e'*velRx;

R(2,2) = 1^2;

H(2,obj.INDS_STATE.VEL) = -[0 1 0]*R_b_e';

measMat(2,:) = [0 0 0 0 pseudoMeasi 4];



end