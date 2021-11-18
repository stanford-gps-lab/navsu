function los_xyzb=findLosXyzb(xyz_usr, xyz_sat, losmask,flagPairs)

%*************************************************************************
%*     Copyright c 2009 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
%FIND_LOS_XYZB calculates the 4D line of sight vectors
%
%LOS_XYZB=FIND_LOS_XYZB(XYZ_USR, XYZ_SAT)
%   Given n_usr user xyz positions and n_sat satellite xyz positions both in
%   ECEF WGS-84 coordinates (X in first column, Y in second column ...) in
%   XYZ_USR and XYZ_SAT respectively, this function returns the n_usr*n_sat
%   line of sight unit vectors from user to satellite augmented by a one at
%   the end.  These LOS vectors may then be used to form the position
%   solution. Optional LOSMASK is a vector of indices (1 to n_usr*n_sat)
%   that specifies which LOS vectors to selectively compute.
%
%   See also: FIND_LOS_ENUB

%2001Mar26 Created by Todd Walter
%2001Apr26 Modified by Wyant Chan   -   Added losmask feature
%2009Nov23 Modified by Todd Walter - Changed sign convention
%2021Oct08 Modified by Fabian Rothmaier - Simplifications and speedup

if nargin < 4
    flagPairs = false;
end

if ~flagPairs
    n_usr = size(xyz_usr, 1);
    n_sat = size(xyz_sat, 1);
    if (nargin==2) || isempty(losmask)
        losmask = true(n_usr * n_sat, 1);
    end
    
    % bring user and satellite positions to the right dimensions
    xyz_usr = repelem(xyz_usr, n_sat, 1);
    xyz_sat = repmat(xyz_sat, n_usr, 1);
    
    % limit to desired lines of sight
    xyz_usr = xyz_usr(losmask, :);
    xyz_sat = xyz_sat(losmask, :);
    
end

% compute xyz lines of sight
los_xyzb = [xyz_sat-xyz_usr ones(size(xyz_usr,1),1)];

%normalize first three columns
mag = sqrt(sum(los_xyzb(:, 1:3).^2, 2));
los_xyzb(:, 1:3) = los_xyzb(:, 1:3) ./ mag;

end