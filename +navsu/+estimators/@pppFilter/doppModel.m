function [predMeas,H,sig] = doppModel(obj,nState,dVel,A,rxDrift,constInd)

% Predicted measurement
predMeas = -dot(dVel,-A)-rxDrift;

% Measurement sensitivity matrix
H = zeros(1,nState);
H(1,obj.INDS_STATE.VEL) = -A;
H(1,obj.INDS_STATE.CLOCK_DRIFT(obj.INDS_STATE.CLOCK_DRIFT_CONSTS == constInd)) = -1;

% Measurement sigma (unweighted by elevation)
sig = obj.PARAMS.sigMeas.dopp;


end