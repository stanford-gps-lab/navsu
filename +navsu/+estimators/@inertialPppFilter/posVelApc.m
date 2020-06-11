function [posOut, velOut] = posVelApc(obj,FLAG_REVERSE)

% Default usage moves the position from the reference position to the
% antenna phase center.  The opposite can also be done.
if nargin < 2
    FLAG_APC_TO_REF = 0;
else
    FLAG_APC_TO_REF = FLAG_REVERSE;
end



% For the standard ppp solution, the reference point is really just the
% antenna phase center anyway.

R_b_e = obj.R_b_e;

rRefApc = obj.PARAMS.ARM_REF_APC;

gyroMeas = obj.lastGyroMeas;
if isempty(gyroMeas)
    gyroMeas = [0 0 0]';
end

Omega_ie = navsu.geo.crossProdMatrix([0,0,navsu.constants.omegaEarth]);

if ~FLAG_APC_TO_REF
    % reference position to antenna phase center
    posOut = obj.pos+R_b_e*rRefApc;
    velOut = obj.vel+R_b_e*navsu.geo.crossProdMatrix(gyroMeas)*rRefApc - ...
        Omega_ie*R_b_e*rRefApc*0;
    
else
    % antenna phase center to reference position
    posOut = obj.pos-R_b_e*rRefApc;
    velOut = obj.vel-R_b_e*navsu.geo.crossProdMatrix(gyroMeas)*rRefApc - ...
        Omega_ie*R_b_e*rRefApc*0;

end


end