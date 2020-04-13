function dcm = euler2dcm123(eulers)

phi   = eulers(1);
theta = eulers(2);
psi   = eulers(3);

% converting from Euler 1-2-3 to direction cosines matrix
% from diebel

% dcm = [cos(phi)*cos(psi)    cos(theta)*sin(psi)    -sin(theta);
%     sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(psi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(theta)*sin(phi);
%     cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(theta)*cos(phi)];

dcm = [cos(theta)*cos(psi)                               cos(theta)*sin(psi)                        -sin(theta);
    sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(psi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(theta)*sin(phi);
    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(theta)*cos(phi)];








end