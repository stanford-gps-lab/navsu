function eulers = dcm2euler123(dcm)

% converting from direction cosines matrix to Euler 1-2-3
% from diebel

phi   = atan2(dcm(2,3),dcm(3,3));
theta = -asin(dcm(1,3));
psi   = atan2(dcm(1,2),dcm(1,1));

eulers = [phi; theta; psi];


end