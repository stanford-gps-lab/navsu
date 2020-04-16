function X_sat_rot = earthRotTravelCorr(traveltime, X_sat)
% earthRotTravelCorr
% DESCRIPTION:
% Returns rotated satellite ECEF coordinates due to Earth rotation during 
% signal travel time
% INPUT:
%   travelTime  - signal travel time
%   X_sat       - satellite's ECEF coordinates
% OUTPUT:
%   X_sat_rot   - rotated satellite's coordinates (ECEF)
%
% See also: 
%
%Written by Kai Borre
%Copyright (c) by Kai Borre
%
% CVS record:
% $Id: e_r_corr.m,v 1.1.1.1.2.6 2006/08/22 13:45:59 dpl Exp $
%==========================================================================

Omegae_dot = 7.292115147e-5;           %  rad/sec

X_sat_rot = zeros(size(X_sat));

for idx = 1:length(traveltime)

%--- Find rotation angle --------------------------------------------------
omegatau   = Omegae_dot * traveltime(idx);

%--- Make a rotation matrix -----------------------------------------------
R3 = [ cos(omegatau)    sin(omegatau)   0;
      -sin(omegatau)    cos(omegatau)   0;
       0                0               1];

%--- Do the rotation ------------------------------------------------------
X_sat_rot(idx,:) = (R3 * X_sat(idx,:)');

end

end