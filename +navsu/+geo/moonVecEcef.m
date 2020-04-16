function [moonVec] = moonVecEcef(jds)
% moonVecEcef
% DESCRIPTION:
%   Produce the lunar position vector in earth centered earth fixed (ECEF)
%   coordinates.  LOS matches to within milliradians to the same output
%   from NASA MICE DE405(or 430?) ephemerides.  This was built from the
%   LunarAzEl function written by Darin Koblick- more information below.
% INPUT:
%   jd      - N-length vector of julian dates
% OUTPUT:
%   moonVec - [Nx3] matrix of lunar positions in ECEF [m]
%
% See also: navsu.time.epochs2jd

%% Revision History:
% Programed by Darin C. Koblick 2/14/2009
%
% Updated on 03/04/2009 to clean up code and add quadrant check to Azimuth
% Thank you Doug W. for your help with the test code to find the quadrant check
% error.
%
% Updated on 04/13/2009 to add Lunar perturbation offsets (this will
% increase the accuracy of calculations)
%
%
% External Source References:
%
% Basics of Positional Astronomy and Ephemerides
% http://jgiesen.de/elevazmoon/basics/index.htm
%
% Computing planetary positions - a tutorial with worked examples
% http://stjarnhimlen.se/comp/tutorial.html
%
%Verified output by comparison with the following source data:
%http://aa.usno.navy.mil/data/docs/AltAz.php

% Code Sequence
%--------------------------------------------------------------------------

nEpochs = length(jds);

moonVec = nan(nEpochs,3);

%Declare Earth Equatorial Radius Measurements in km
EarthRadEq = 6378.1370/1000;

for idx = 1:nEpochs
    jd = jds(idx);
    %Find the Day Number
    d = jd - 2451543.5;
    
    %Keplerian Elements of the Moon
    %This will also account for the Sun's perturbation
    N = 125.1228-0.0529538083.*d; %    (Long asc. node deg)
    i = 5.1454; %                      (Inclination deg)
    w = 318.0634 + 0.1643573223.*d; %  (Arg. of perigee deg)
    a =  60.2666;%                     (Mean distance (Earth's Equitorial Radii)
    e = 0.054900;%                     (Eccentricity)
    M = mod(115.3654+13.0649929509.*d,360);%    (Mean anomaly deg)
    
    LMoon =  mod(N + w + M,360);                 %(Moon's mean longitude deg)
    FMoon =  mod(LMoon - N,360);                 %(Moon's argument of latitude)
    
    L = w+M;
    
    %Keplerian Elements of the Sun
    wSun = mod(282.9404 + 4.70935E-5.*d,360);    % (longitude of perihelion)
    MSun = mod(356.0470 + 0.9856002585.*d,360);  % (Sun mean anomaly)
    LSun = mod(wSun + MSun,360);                 % (Sun's mean longitude)
    
    DMoon =  LMoon - LSun;                     % (Moon's mean elongation)
    
    
    %Calculate Lunar perturbations in Longitude
    LunarPLon = [ -1.274.*sin((M - 2.*DMoon).*(pi/180)); ...
        .658.*sin(2.*DMoon.*(pi/180)); ...
        -0.186.*sin(MSun.*(pi/180)); ...
        -0.059.*sin((2.*M-2.*DMoon).*(pi/180)); ...
        -0.057.*sin((M-2.*DMoon + MSun).*(pi/180)); ...
        .053.*sin((M+2.*DMoon).*(pi/180)); ...
        .046.*sin((2.*DMoon-MSun).*(pi/180)); ...
        .041.*sin((M-MSun).*(pi/180)); ...
        -0.035.*sin(DMoon.*(pi/180)); ...
        -0.031.*sin((M+MSun).*(pi/180)); ...
        -0.015.*sin((2.*FMoon-2.*DMoon).*(pi/180)); ...
        .011.*sin((M-4.*DMoon).*(pi/180))];
    
    %Calculate Lunar perturbations in Latitude
    LunarPLat = [ -0.173.*sin((FMoon-2.*DMoon).*(pi/180)); ...
        -0.055.*sin((M-FMoon-2.*DMoon).*(pi/180)); ...
        -0.046.*sin((M+FMoon-2.*DMoon).*(pi/180)); ...
        +0.033.*sin((FMoon+2.*DMoon).*(pi/180)); ...
        +0.017.*sin((2.*M+FMoon).*(pi/180))];
    
    %Calculate perturbations in Distance
    LunarPDist = [ -0.58*cos((M-2.*DMoon).*(pi/180)); ...
        -0.46.*cos(2.*DMoon.*(pi/180))];
    
    % Compute E, the eccentric anomaly
    
    %E0 is the eccentric anomaly approximation estimate
    %(this will initially have a relativly high error)
    E0 = M+(180./pi).*e.*sin(M.*(pi/180)).*(1+e.*cos(M.*(pi/180)));
    
    %Compute E1 and set it to E0 until the E1 == E0
    E1 = E0-(E0-(180/pi).*e.*sin(E0.*(pi/180))-M)./(1-e*cos(E0.*(pi/180)));
    while E1-E0 > .000005
        E0 = E1;
        E1 = E0-(E0-(180/pi).*e.*sin(E0.*(pi/180))-M)./(1-e*cos(E0.*(pi/180)));
    end
    E = E1;
    
    %Compute rectangular coordinates (x,y) in the plane of the lunar orbit
    x = a.*(cos(E.*(pi/180))-e);
    y = a.*sqrt(1-e.*e).*sin(E.*(pi/180));
    
    %convert this to distance and true anomaly
    r = sqrt(x.*x + y.*y);
    v = atan2(y.*(pi/180),x.*(pi/180)).*(180/pi);
    
    %Compute moon's position in ecliptic coordinates
    xeclip = r.*(cos(N.*(pi/180)).*cos((v+w).*(pi/180))-sin(N.*(pi/180)).*sin((v+w).*(pi/180)).*cos(i.*(pi/180)));
    yeclip = r.*(sin(N.*(pi/180)).*cos((v+w).*(pi/180))+cos(N.*(pi/180))*sin(((v+w).*(pi/180)))*cos(i.*(pi/180)));
    zeclip = r.*sin((v+w).*(pi/180)).*sin(i.*(pi/180));
    
    %Add the calculated lunar perturbation terms to increase model fidelity
    [eLon, eLat, eDist] = cart2sph(xeclip,yeclip,zeclip);
    [xeclip, yeclip, zeclip] = sph2cart(eLon + sum(LunarPLon).*(pi/180), ...
        eLat + sum(LunarPLat).*(pi/180), ...
        eDist + sum(LunarPDist));
    
    %convert the latitude and longitude to right ascension RA and declination
    %delta
    T = (jd-2451545.0)/36525.0;
    
    %Generate a rotation matrix for ecliptic to equitorial
    %RotM=rotm_coo('E',jd);
    %See rotm_coo.m for obl and rotational matrix transformation
    Obl = 23.439291 - 0.0130042.*T - 0.00000016.*T.*T + 0.000000504.*T.*T.*T;
    Obl = Obl.*(pi/180);
    
    RotM = [1 0 0; 0 cos(Obl) sin(Obl); 0 -sin(Obl) cos(Obl)]';
    
    %Apply the rotational matrix to the ecliptic rectangular coordinates
    %Also, convert units to km instead of earth equatorial radii
    sol = RotM*[xeclip yeclip zeclip]'.*EarthRadEq*1000*1000;
    
    xequat = sol(1);
    yequat = sol(2);
    zequat = sol(3);
    
    %convert equatorial rectangular coordinates to RA and Decl:
    Lon = 0;
    
    r = sqrt(xequat.^2 + yequat.^2 + zequat.^2);
    RA = atan2(yequat,xequat).*(180/pi);
    delta = asin(zequat./r).*(180/pi);
    
    %Following the RA DEC to Az Alt conversion sequence explained here:
    %http://www.stargazing.net/kepler/altaz.html
    
    %Find the J2000 value
    J2000 = jd - 2451545.0;
    hourvec = datevec(navsu.time.epochs2datenum(navsu.time.jd2epochs(jd)));
    UTH = hourvec(:,4) + hourvec(:,5)/60 + hourvec(:,6)/3600;
    
    %Calculate local siderial time
    GMST0=mod(L+180,360)./15;
    
    LST = mod(100.46+0.985647.*J2000+Lon+15*UTH,360);
    
    theta = -LST*pi/180;
    posEcef = [xequat.*cos(theta)-yequat.*sin(theta)    xequat.*sin(theta)+yequat.*cos(theta) zequat];
    
    moonVec(idx,:) = posEcef;
   
end





