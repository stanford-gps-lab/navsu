%
function xcorsta = step2diu(xsta,fhr,t)
%+
%  - - - - - - - - - - -
%   S T E P 2 D I U
%  - - - - - - - - - - -
%
%  This routine is part of the International Earth Rotation and
%  Reference Systems Service (IERS) Conventions software collection.
%
%  This subroutine gives the in-phase and out-of-phase corrections
%  induced by mantle anelasticity in the diurnal band.
%
%  In general, Class 1, 2, and 3 models represent physical effects that
%  act on geodetic parameters while canonical models provide lower-level
%  representations or basic computations that are used by Class 1, 2, or
%  3 models.
%
%  Status: Class 1
%
%     Class 1 models are those recommended to be used a priori in the
%     reduction of raw space geodetic data in order to determine
%     geodetic parameter estimates.
%     Class 2 models are those that eliminate an observational
%     singularity and are purely conventional in nature.
%     Class 3 models are those that are not required as either Class
%     1 or 2.
%     Canonical models are accepted as is and cannot be classified as a
%     Class 1, 2, or 3 model.
%
%  Given:
%     XSTA          d(3)   Geocentric position of the IGS station (Note 1)
%     FHR           d      Fractional hours in the day (Note 2)
%     T             d      Centuries since J2000
%
%  Returned:
%     XCORSTA       d(3)   In phase and out of phase station corrections
%                          for diurnal band (Note 4)
%
%  Notes:
%
%  1) The IGS station is in ITRF co-rotating frame.  All coordinates are
%     expressed in meters.
%
%  2) The fractional hours in the day is computed as the hour + minutes/60.0
%     + sec/3600.0.  The unit is expressed in Universal Time (UT).
%
%  4) All coordinates are expressed in meters.
%
%  Test case:
%     given input: XSTA(1) = 4075578.385D0 meters
%                  XSTA(2) =  931852.890D0 meters
%                  XSTA(3) = 4801570.154D0 meters
%                  FHR     = 0.00D0 hours
%                  T       = 0.1059411362080767D0 Julian centuries
%
%     expected output:  XCORSTA(1) = 0.4193085327321284701D-02 meters
%                       XCORSTA(2) = 0.1456681241014607395D-02 meters
%                       XCORSTA(3) = 0.5123366597450316508D-02 meters
%
%  References:
%
%     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
%     displacements,' J. Geophys. Res., 102(B9), pp. 20,469-20,477
%
%     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
%     IERS Technical Note No. 36, BKG (2010)
%
%  Revisions:
%  1996 March    23 V. Dehant      Original code
%  2009 July     31 B.E. Stetzler  Initial standardization of code
%  2009 August   06 B.E. Stetzler  Provided a test case
%  2009 August   06 B.E. Stetzler  Capitalized all variables for
%                                  Fortran 77 compatibility
%  2010 October  20 B.E. Stetzler  Input T corrected to be number of
%                                  centuries since J2000
%-----------------------------------------------------------------------


d2pi=6.283185307179586476925287d0;

[ datdi([1:9],[1:31])]=reshape([-3d0,0d0,2d0,0d0,0d0,-0.01d0,0d0,0d0,0d0,-3d0,2d0,0d0,0d0,0d0,-0.01d0,0d0,0d0,0d0,-2d0,0d0,1d0,-1d0,0d0,-0.02d0,0d0,0d0,0d0,-2d0,0d0,1d0,0d0,0d0,-0.08d0,0d0,-0.01d0,0.01d0,-2d0,2d0,-1d0,0d0,0d0,-0.02d0,0d0,0d0,0d0,-1d0,0d0,0d0,-1d0,0d0,-0.10d0,0d0,0d0,0d0,-1d0,0d0,0d0,0d0,0d0,-0.51d0,0d0,-0.02d0,0.03d0,-1d0,2d0,0d0,0d0,0d0,0.01d0,0d0,0d0,0d0,0d0,-2d0,1d0,0d0,0d0,0.01d0,0d0,0d0,0d0,0d0,0d0,-1d0,0d0,0d0,0.02d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0.06d0,0d0,0d0,0d0,0d0,0d0,1d0,1d0,0d0,0.01d0,0d0,0d0,0d0,0d0,2d0,-1d0,0d0,0d0,0.01d0,0d0,0d0,0d0,1d0,-3d0,0d0,0d0,1d0,-0.06d0,0d0,0d0,0d0,1d0,-2d0,0d0,-1d0,0d0,0.01d0,0d0,0d0,0d0,1d0,-2d0,0d0,0d0,0d0,-1.23d0,-0.07d0,0.06d0,0.01d0,1d0,-1d0,0d0,0d0,-1d0,0.02d0,0d0,0d0,0d0,1d0,-1d0,0d0,0d0,1d0,0.04d0,0d0,0d0,0d0,1d0,0d0,0d0,-1d0,0d0,-0.22d0,0.01d0,0.01d0,0d0,1d0,0d0,0d0,0d0,0d0,12.00d0,-0.80d0,-0.67d0,-0.03d0,1d0,0d0,0d0,1d0,0d0,1.73d0,-0.12d0,-0.10d0,0d0,1d0,0d0,0d0,2d0,0d0,-0.04d0,0d0,0d0,0d0,1d0,1d0,0d0,0d0,-1d0,-0.50d0,-0.01d0,0.03d0,0d0,1d0,1d0,0d0,0d0,1d0,0.01d0,0d0,0d0,0d0,0d0,1d0,0d0,1d0,-1d0,-0.01d0,0d0,0d0,0d0,1d0,2d0,-2d0,0d0,0d0,-0.01d0,0d0,0d0,0d0,1d0,2d0,0d0,0d0,0d0,-0.11d0,0.01d0,0.01d0,0d0,2d0,-2d0,1d0,0d0,0d0,-0.01d0,0d0,0d0,0d0,2d0,0d0,-1d0,0d0,0d0,-0.02d0,0d0,0d0,0d0,3d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,3d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0,0d0],9,31);

deg2rad = d2pi./360d0;

%  Compute the phase angles in degrees.
s = 218.31664563d0 +(481267.88194d0+(-0.0014663889d0+(0.00000185139d0).*t).*t).*t;

tau = fhr.*15d0 + 280.4606184d0 +(36000.7700536d0+(0.00038793d0+(-0.0000000258d0).*t).*t).*t +(-s);

pr =(1.396971278d0+(0.000308889d0+(0.000000021d0+(0.000000007d0).*t).*t).*t).*t;

s = s + pr;

h = 280.46645d0 +(36000.7697489d0+(0.00030322222d0+(0.000000020d0+(-0.00000000654d0).*t).*t).*t).*t;

p = 83.35324312d0 +(4069.01363525d0+(-0.01032172222d0+(-0.0000124991d0+(0.00000005263d0).*t).*t).*t).*t;

zns = 234.95544499d0 +(1934.13626197d0+(-0.00207561111d0+(-0.00000213944d0+(0.00000001650d0).*t).*t).*t).*t;

ps = 282.93734098d0 +(1.71945766667d0+(0.00045688889d0+(-0.00000001778d0+(-0.00000000334d0).*t).*t).*t).*t;

% Reduce angles to between the range 0 and 360.
s = rem(s,360d0);
tau = rem(tau,360d0);
h = rem(h,360d0);
p = rem(p,360d0);
zns = rem(zns,360d0);
ps = rem(ps,360d0);

rsta = sqrt(xsta(1).^2+xsta(2).^2+xsta(3).^2);
sinphi = xsta(3)./rsta;
cosphi = sqrt(xsta(1).^2+xsta(2).^2)./rsta;

cosla = xsta(1)./cosphi./rsta;
sinla = xsta(2)./cosphi./rsta;
zla = atan2(xsta(2),xsta(1));

for i = 1 : 3
    % Initialize.
    xcorsta(i) = 0d0;
end; i = fix(3+1);
for j = 1 : 31
    % Convert from degrees to radians.
    thetaf =(tau+datdi(1,j).*s+datdi(2,j).*h+datdi(3,j).*p+datdi(4,j).*zns+datdi(5,j).*ps).*deg2rad;
    
    dr = datdi(6,j).*2d0.*sinphi.*cosphi.*sin(thetaf+zla) + datdi(7,j).*2d0.*sinphi.*cosphi.*cos(thetaf+zla);
    
    dn = datdi(8,j).*(cosphi.^2-sinphi.^2).*sin(thetaf+zla)+ datdi(9,j).*(cosphi.^2-sinphi.^2).*cos(thetaf+zla);
    %     Modified 20 June 2007
    
    de = datdi(8,j).*sinphi.*cos(thetaf+zla) - datdi(9,j).*sinphi.*sin(thetaf+zla);
    
    xcorsta(1) = xcorsta(1) + dr.*cosla.*cosphi - de.*sinla -dn.*sinphi.*cosla;
    xcorsta(2) = xcorsta(2) + dr.*sinla.*cosphi + de.*cosla -dn.*sinphi.*sinla;
    xcorsta(3) = xcorsta(3) + dr.*sinphi + dn.*cosphi;
end; j = fix(31+1);

for i = 1 : 3
    xcorsta(i) = xcorsta(i)./1000d0;
end; i = fix(3+1);

%  Finished.

%+----------------------------------------------------------------------
%
%  Copyright (C) 2008
%  IERS Conventions Center
%
%  ==================================
%  IERS Conventions Software License
%  ==================================
%
%  NOTICE TO USER:
%
%  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
%  WHICH APPLY TO ITS USE.
%
%  1. The Software is provided by the IERS Conventions Center ('the
%     Center').
%
%  2. Permission is granted to anyone to use the Software for any
%     purpose, including commercial applications, free of charge,
%     subject to the conditions and restrictions listed below.
%
%  3. You (the user) may adapt the Software and its algorithms for your
%     own purposes and you may distribute the resulting 'derived work'
%     to others, provided that the derived work complies with the
%     following requirements:
%
%     a) Your work shall be clearly identified so that it cannot be
%        mistaken for IERS Conventions software and that it has been
%        neither distributed by nor endorsed by the Center.
%
%     b) Your work (including source code) must contain descriptions of
%        how the derived work is based upon and/or differs from the
%        original Software.
%
%     c) The name(s) of all modified routine(s) that you distribute
%        shall be changed.
%
%     d) The origin of the IERS Conventions components of your derived
%        work must not be misrepresented; you must not claim that you
%        wrote the original Software.
%
%     e) The source code must be included for all routine(s) that you
%        distribute.  This notice must be reproduced intact in any
%        source distribution.
%
%  4. In any published work produced by the user and which includes
%     results achieved by using the Software, you shall acknowledge
%     that the Software was used in obtaining those results.
%
%  5. The Software is provided to the user 'as is' and the Center makes
%     no warranty as to its use or performance.   The Center does not
%     and cannot warrant the performance or results which the user may
%     obtain by using the Software.  The Center makes no warranties,
%     express or implied, as to non-infringement of third party rights,
%     merchantability, or fitness for any particular purpose.  In no
%     event will the Center be liable to the user for any consequential,
%     incidental, or special damages, including any lost profits or lost
%     savings, even if a Center representative has been advised of such
%     damages, or for any claim by any third party.
%
%  Correspondence concerning IERS Conventions software should be
%  addressed as follows:
%
%                     Gerard Petit
%     Internet email: gpetit[at]bipm.org
%     Postal address: IERS Conventions Center
%                     Time, frequency and gravimetry section, BIPM
%                     Pavillon de Breteuil
%                     92312 Sevres  FRANCE
%
%     or
%
%                     Brian Luzum
%     Internet email: brian.luzum[at]usno.navy.mil
%     Postal address: IERS Conventions Center
%                     Earth Orientation Department
%                     3450 Massachusetts Ave, NW
%                     Washington, DC 20392
%
%
%-----------------------------------------------------------------------
end

