function [posApc, velApc] = posVelApc(obj)


% For the standard ppp solution, the reference point is really just the
% antenna phase center anyway. 

% TO DO: MAKE THIS MORE FLEXIBLE- SEE HOW IT IS DONE IN THE
% INERTIALPPPFILTER AND MAYBE GO FROM THERE?

posApc = obj.pos;
velApc = obj.vel;



end