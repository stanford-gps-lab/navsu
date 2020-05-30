function [posApc, velApc] = posVelApc(obj)


% For the standard ppp solution, the reference point is really just the
% antenna phase center anyway. 



R_b_e = obj.R_b_e;

rRefApc = obj.PARAMS.IMU_ARM;

gyroMeas = obj.lastGyroMeas;
if isempty(gyroMeas)
   gyroMeas = [0 0 0]'; 
end

Omega_ie = navsu.geo.crossProdMatrix([0,0,navsu.constants.omegaEarth]);

% 

posApc = obj.pos;
velApc = obj.vel;

posApc = obj.pos+R_b_e*rRefApc;
velApc = obj.vel+R_b_e*navsu.geo.crossProdMatrix(gyroMeas)*rRefApc - ...
    Omega_ie*R_b_e*rRefApc*0;


end