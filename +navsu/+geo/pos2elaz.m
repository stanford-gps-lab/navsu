function [el,az] = pos2elaz(usrPos,satPos)

llh = navsu.geo.xyz2llh(usrPos);

RxyzEnu = navsu.geo.findxyz2enu(llh(1)*pi/180,llh(2)*pi/180);
usr_ehat = RxyzEnu(:,1);
usr_nhat = RxyzEnu(:,2);
usr_uhat = RxyzEnu(:,3);


losxyzb = findLosXyzb(usrPos,satPos);

los_enub = calcLosEnub(losxyzb,usr_ehat',usr_nhat',usr_uhat');
[el, az] = findElaz(los_enub);


end


function [el, az]=findElaz(los_enub)

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
%FIND_ELAZ calculates the elevation and azimuth angles from the ENU LOS vectors
%
%[EL AZ]=FIND_ELAZ(LOS_ENU);
%   Given line of sight vectors in the East, North, Up frame this function
%   calculates the elevation and azimuth angles in radians.
%
%   See also: FIND_LOS_XYZB

%2001Mar26 Created by Todd Walter


el=asin(-los_enub(:,3));
az=atan2(-los_enub(:,1),-los_enub(:,2));
end

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
%   line of sight unit vectors augmented by a one at the end.  These LOS
%   vectors may then be used to form the position solution.  Optional LOSMASK
%   is a vector of indices (1 to n_usr*n_sat) that specifies which LOS vectors
%   to selectively compute.
%
%   See also: FIND_LOS_ENUB

%2001Mar26 Created by Todd Walter
%2001Apr26 Modified by Wyant Chan   -   Added losmask feature
%2009Nov23 Modified by Todd Walter - Changed sign convention

if nargin < 4
    flagPairs = 0;
end

if ~flagPairs
    [n_usr tmp]=size(xyz_usr);
    [n_sat tmp]=size(xyz_sat);
    n_los=n_usr*n_sat;
    if (nargin==2)
        losmask = [1:n_los]';
    end
    n_mask = size(losmask,1);
    
    %initialize 4th column of the line of sight vector
    los_xyzb = ones(n_mask,4);
    
    %build the line of sight vector
    [t1 t2]=meshgrid(xyz_usr(:,1),xyz_sat(:,1));
    t1 = reshape(t1,n_los,1);
    t2 = reshape(t2,n_los,1);
    los_xyzb(:,1) = t2(losmask) - t1(losmask);
    
    [t1 t2]=meshgrid(xyz_usr(:,2),xyz_sat(:,2));
    t1 = reshape(t1,n_los,1);
    t2 = reshape(t2,n_los,1);
    los_xyzb(:,2) = t2(losmask) - t1(losmask);
    
    [t1 t2]=meshgrid(xyz_usr(:,3),xyz_sat(:,3));
    t1 = reshape(t1,n_los,1);
    t2 = reshape(t2,n_los,1);
    los_xyzb(:,3) = t2(losmask) - t1(losmask);
    
    %normalize first three columns
    mag=sqrt(sum(los_xyzb(:,1:3)'.^2))';
    los_xyzb(:,1)=los_xyzb(:,1)./mag;
    los_xyzb(:,2)=los_xyzb(:,2)./mag;
    los_xyzb(:,3)=los_xyzb(:,3)./mag;
    
else
    
    los_xyzb = [xyz_sat-xyz_usr ones(size(xyz_usr,1),1)];
    
    %normalize first three columns
    mag=sqrt(sum(los_xyzb(:,1:3)'.^2))';
    los_xyzb(:,1)=los_xyzb(:,1)./mag;
    los_xyzb(:,2)=los_xyzb(:,2)./mag;
    los_xyzb(:,3)=los_xyzb(:,3)./mag;
    
end
end

function los_enub = calcLosEnub(los_xyzb, e_hat, n_hat, u_hat)

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
%CALC_LOS_ENUB converts 4D line of sight vectors from XYZ to East North Up
%
%LOS_ENUB=CALC_LOS_XYZB(LOS_XYZB, E_HAT, N_HAT, U_HAT);
%   Given n_los line of sight vectors in ECEF WGS-84 coordinates (X in first
%   column, Y in second column, Z in the third column and 1 in the fourth) in
%   LOS_XYZB and n_los east, north, and up unit vectors (at the user location)
%   in E_HAT, N_HAT, and U_HAT respectively this function returns the n_los
%   line of sight unit vectors augmented by a one at the end in the East, North
%   and Up frame.  These LOS vectors may then be used to form the position
%   solution.
%
%   See also: FIND_LOS_XYZB FIND_LOS_ENUB

%2001Mar26 Created by Todd Walter
%2009Nov23 Modified by Todd Walter - Changed sign convention

% Check if e,n,u_hat are correct size and repmat if not
if size(e_hat,1) == 1 && size(los_xyzb,1) ~= 1
    e_hat = repmat(e_hat,size(los_xyzb,1),1);
    n_hat = repmat(n_hat,size(los_xyzb,1),1);
    u_hat = repmat(u_hat,size(los_xyzb,1),1);
end

%initialize 4th column of the line of sight vector
los_enub(:,4)=los_xyzb(:,4);

%dot the east unit vector with the los vector to determine -cos(elev)*sin(azim)
los_enub(:,1)=sum((-e_hat.*los_xyzb(:,1:3))')';

%dot the north unit vector with the los vector to determine -cos(elev)cos(azim)
los_enub(:,2)=sum((-n_hat.*los_xyzb(:,1:3))')';

%dot the up unit vector with the los vector to determine -sin(elevation)
los_enub(:,3)=sum((-u_hat.*los_xyzb(:,1:3))')';

end




