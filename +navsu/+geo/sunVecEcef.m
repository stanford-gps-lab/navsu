function sunVec = sunVecEcef(jd)
% sunVecEcef
% DESCRIPTION:
%   Produce the solar position vector in earth centered earth fixed (ECEF)
%   coordinates.  LOS matches to within milliradians to the same output
%   from NASA MICE DE405(or 430?) ephemerides.  This was built from the
%   SolarAzEl function written by Darin Koblick- more information below. 
% INPUT:
%   jd     - N-length vector of julian dates
% OUTPUT:
%   sunVec - [Nx3] matrix of solar position vectors in ECEF [m]
%
% See also: navsu.time.epochs2jd

%% Revision History:
% Programed by Darin C. Koblick 2/17/2009
%
%              Darin C. Koblick 4/16/2013 Vectorized for Speed
%                                         Allow for MATLAB Datevec input in
%                                         addition to a UTC string.
%                                         Cleaned up comments and code to
%                                         avoid warnings in MATLAB editor.
%
%--------------------------------------------------------------------------
% External Function Call Sequence:
%[Az El] = SolarAzEl('1991/05/19 13:00:00',50,10,0)
%% Output Description:
%Az                     [N x 1]         Azimuth location of the sun (deg)
%El                     [N x 1]         Elevation location of the sun (deg)
%
%
%% Source References:
%Solar Position obtained from:
%http://stjarnhimlen.se/comp/tutorial.html#5
%% Begin Code Sequence
% Ensure that input is a column vector 
jd = jd(:);

d = jd-2451543.5;

% Keplerian Elements for the Sun (geocentric)
w = 282.9404+4.70935e-5*d; %    (longitude of perihelion degrees)
a = 1.495978707e11;              %(mean distance, m)
e = 0.016709-1.151e-9.*d;%       (eccentricity)
M = mod(356.0470+0.9856002585.*d,360);%   (mean anomaly degrees)
L = w + M;                     %(Sun's mean longitude degrees)
oblecl = 23.4393-3.563e-7.*d;  %(Sun's obliquity of the ecliptic)

%auxiliary angle
E = M+(180/pi).*e.*sin(M.*(pi/180)).*(1+e.*cos(M.*(pi/180)));

%rectangular coordinates in the plane of the ecliptic (x axis toward
%perhilion)
x = a*(cos(E.*(pi/180))-e);
y = a*sin(E.*(pi/180)).*sqrt(1-e.^2);

%find the distance and true anomaly
r = sqrt(x.^2 + y.^2);
v = atan2(y,x).*(180/pi);

%find the longitude of the sun
lon = v + w;

%compute the ecliptic rectangular coordinates
xeclip = r.*cos(lon.*(pi/180));
yeclip = r.*sin(lon.*(pi/180));
zeclip = 0.0;

%rotate these coordinates to equitorial rectangular coordinates
xequat = xeclip;
yequat = yeclip.*cos(oblecl.*(pi/180))+zeclip*sin(oblecl.*(pi/180));
zequat = yeclip.*sin(23.4406.*(pi/180))+zeclip*cos(oblecl.*(pi/180));

%convert equatorial rectangular coordinates to RA and Decl:
Lon = 0;

r = sqrt(xequat.^2 + yequat.^2 + zequat.^2); 
RA = atan2(yequat,xequat).*(180/pi);
delta = asin(zequat./r).*(180/pi);

%Following the RA DEC to Az Alt conversion sequence explained here:
%http://www.stargazing.net/kepler/altaz.html

%Find the J2000 value
%J2000 = jd - 2451545.0;
hourvec = datevec(navsu.time.epochs2datenum(navsu.time.jd2epochs(jd)));
UTH = hourvec(:,4) + hourvec(:,5)/60 + hourvec(:,6)/3600;

%Calculate local siderial time
GMST0=mod(L+180,360)./15;
SIDTIME = GMST0 + UTH + Lon./15;

theta   = -SIDTIME.*15*pi/180;
posEcef = [xequat.*cos(theta)-yequat.*sin(theta)    xequat.*sin(theta)+yequat.*cos(theta) zequat];

sunVec = posEcef;
end

















