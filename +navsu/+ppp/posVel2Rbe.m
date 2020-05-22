function R_b_e = posVel2Rbe(pos,vel)


% Given a position and velocity in ECEF, output a rotation matrix from body
% frame to ECEF. 

% x is backwards
% y is rightward
% z is up

xi  = -vel./norm(vel);
z0i = pos./norm(pos);
yi  = -cross(xi,z0i);
zi = cross(xi,yi);

R_b_e = [xi yi zi];



end