function inertialNavEquationsEcef(obj,epoch,accMeasi,gyroMeasi)


oldPos = obj.pos;
oldVel = obj.vel;
old_R_b_e = obj.R_b_e;

dt = min([epoch-obj.epochLastInertialUpdate epoch-obj.epochLastGnssUpdate]);
% dt = 0.01;
obj.epochLastInertialUpdate = epoch;

% parameters
omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s)

% Begins

% ATTITUDE UPDATE
% From (2.145) determine the Earth rotation over the update interval
% C_Earth = C_e_i' * old_C_e_i
alpha_ie = omega_ie * dt;
C_Earth = [cos(alpha_ie), sin(alpha_ie), 0;...
          -sin(alpha_ie), cos(alpha_ie), 0;...
                       0,             0,  1];

% Calculate attitude increment, magnitude, and skew-symmetric matrix

alpha_ib_b = gyroMeasi * dt;
mag_alpha = sqrt(alpha_ib_b' * alpha_ib_b);
Alpha_ib_b = crossProdMatrix(alpha_ib_b);  

% Obtain coordinate transformation matrix from the new attitude w.r.t. an
% inertial frame to the old using Rodrigues' formula, (5.73)
if mag_alpha>1.E-8
    C_new_old = eye(3) + sin(mag_alpha) / mag_alpha * Alpha_ib_b +...
        (1 - cos(mag_alpha)) / mag_alpha^2 * Alpha_ib_b * Alpha_ib_b;
else
    C_new_old = eye(3) + Alpha_ib_b;
end %if mag_alpha    

% Update attitude using (5.75)
R_b_e = C_Earth * old_R_b_e * C_new_old;

% SPECIFIC FORCE FRAME TRANSFORMATION
% Calculate the average body-to-ECEF-frame coordinate transformation
% matrix over the update interval using (5.84) and (5.85)
if mag_alpha>1.E-8
    ave_C_b_e = old_R_b_e * (eye(3) + (1 - cos(mag_alpha)) / mag_alpha^2 ...
        * Alpha_ib_b + (1 - sin(mag_alpha) / mag_alpha) / mag_alpha^2 ...
        * Alpha_ib_b * Alpha_ib_b) - 0.5 * crossProdMatrix([0;0;alpha_ie])...
        * old_R_b_e;
else
     ave_C_b_e = old_R_b_e - 0.5 * crossProdMatrix([0;0;alpha_ie]) *...
         old_R_b_e;
end %if mag_alpha     

% Transform specific force to ECEF-frame resolving axes using (5.85)
f_ib_e = ave_C_b_e * accMeasi;

% UPDATE VELOCITY
% From (5.36),
vel = oldVel + dt * (f_ib_e + Gravity_ECEF(oldPos) -...
    2 * crossProdMatrix([0;0;omega_ie]) * oldVel);
% UPDATE CARTESIAN POSITION
% From (5.38),
pos = oldPos + (vel + oldVel) * 0.5 * dt; 



%% Update object
obj.pos = pos;
obj.vel = vel;
obj.R_b_e = R_b_e;

obj.lastAccMeas = accMeasi;
obj.lastGyroMeas = gyroMeasi;

end






















