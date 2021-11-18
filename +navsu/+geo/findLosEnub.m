function los_enub=findLosEnub(los_xyzb, usr_ehat, usr_nhat, usr_uhat, losmask)

%*************************************************************************
%*     Copyright c 2001 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
%FIND_LOS_ENUB converts 4D line of sight vectors from XYZ to East North Up
%
%LOS_ENUB=FIND_LOS_XYZB(LOS_XYZB, USR_EHAT, USR_NHAT, USR_UHAT);
%   Given n_los line of sight vectors in ECEF WGS-84 coordinates (X in first
%   column, Y in second column, Z in the third column and 1 in the fourth) in
%   LOS_XYZB and n_usr east, north, and up unit vectors (at the user location)
%   in E_HAT, N_HAT, and U_HAT respectively this function returns the n_los
%   line of sight unit vectors augmented by a one at the end in the East, North
%   and Up frame.  These LOS vectors may then be used to form the position
%   solution.  Optional LOSMASK is a vector of indices (1 to n_usr*n_sat) 
%   that specifies which LOS vectors to selectively compute.
%  
%   See also: FIND_LOS_XYZB CALC_LOS_ENUB

%2001Mar26 Created by Todd Walter
%2001Apr26 Modified by Wyant Chan   -   Added losmask feature
%2021Oct08 Modified by Fabian Rothmaier - simplified the code

n_los = size(los_xyzb, 1);
n_usr = size(usr_ehat, 1);
n_sat = n_los/n_usr;
if nargin == 4
    losmask = true(n_los, 1);
end
sat_idx = 1:n_sat;

%expand the user east unit vector to match the lines of sight
[t1, ~] = meshgrid(usr_ehat, sat_idx);
e_hat = reshape(t1, n_los, 3);
% this is essentially repelem(usr_ehat, n_sat, 1) that works before Matlab
% R2015a

%expand the user north unit vector to match the lines of sight
[t1, ~] = meshgrid(usr_nhat, sat_idx);
n_hat = reshape(t1, n_los, 3);

%expand the user up unit vector to match the lines of sight
[t1, ~] = meshgrid(usr_uhat, sat_idx);
u_hat = reshape(t1, n_los, 3);

%calculate the LOS vectors in the ENU frame
los_enub = navsu.geo.calcLosEnub(los_xyzb(losmask, :), ...
                                 e_hat(losmask, :), ...
                                 n_hat(losmask, :), ...
                                 u_hat(losmask, :));
    
end




