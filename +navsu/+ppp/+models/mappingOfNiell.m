function [Mw,Md] = mappingOfNiell(el, h, lat,doy)

% Source: http://www.navipedia.net/index.php/Mapping_of_Niell
% this also just comes from the original paper: 
% Global mapping functions for the atmosphere delay at radio wavelengths
% (Niell 1996)

% Two columns- hydrostatic and wet

latTable = [15 30 45 60 75];
elRad = el*pi/180;
%% Wet mapping function
% Interpolate a, b, c

wetTable = [...
% lat:  15           30            45          60            75
    5.8021897e-4 5.6794847e-4 5.8118017e-4 5.9727542e-4 6.1641693e-4; ... % a
    1.4275268e-3 1.5138625e-3 1.4572752e-3 1.5007428e-3 1.7599082e-3; ... % b
    4.3472961e-2 4.6729510e-2 4.3908931e-2 4.4626982e-2 5.4736038e-2];    % c

abc = zeros(length(lat),3);
for idx = 1:3
    abc(:,idx)  = interp1(latTable,wetTable(idx,:),lat);
   
    abc(lat <= latTable(1),idx) = wetTable(idx,1);
    abc(lat >= latTable(end),idx) = wetTable(idx,end);
end
a = abc(:,1);
b = abc(:,2);
c = abc(:,3);

Mw = (1+a./(1+b./(1+c)))./(sin(elRad)+a./(sin(elRad)+b./(sin(elRad)+c)));


%% Dry mapping function
dryTable = [ ...
% lat: 
    1.2769934e-3 1.2683230e-3 1.2465397e-3 1.2196049e-3 1.2045996e-3;
    2.9153695e-3 2.9152299e-3 2.9288445e-3 2.9022565e-3 2.9024912e-3;
    62.610505e-3 62.837393e-3 63.721774e-3 63.824265e-3 64.258455e-3;
    0            1.2709626e-5 2.6523662e-5 3.4000452e-5 4.1202191e-5;
    0            2.1414979e-5 3.0160779e-5 7.2562722e-5 11.723375e-5;
    0            9.0128400e-5 4.3497037e-5 84.795348e-5 170.37206e-5];

heightCorrAbc = [2.53e-5 5.49e-3 1.14e-3];

abc2 = zeros(length(lat),6);
for idx = 1:6
    abc2(:,idx)  = interp1(latTable,dryTable(idx,:),lat);
   
    % Areas lower than 15 deg or greater than 75 are just set to the first
    % or last value, respectively
    abc2(lat <= latTable(1),idx)   = dryTable(idx,1);
    abc2(lat >= latTable(end),idx) = dryTable(idx,end);
end

abcDry = zeros(length(lat),3);
aDry = abc2(:,1)-abc2(:,4).*cos(2*pi*(doy-28)/365.25);
bDry = abc2(:,2)-abc2(:,5).*cos(2*pi*(doy-28)/365.25);
cDry = abc2(:,3)-abc2(:,6).*cos(2*pi*(doy-28)/365.25);

Md0 = (1+aDry./(1+bDry./(1+cDry)))./(sin(elRad)+aDry./(sin(elRad)+bDry./(sin(elRad)+cDry)));

aHt = heightCorrAbc(1);
bHt = heightCorrAbc(2);
cHt = heightCorrAbc(3);

mht = (1+aHt./(1+bHt./(1+cHt)))./(sin(elRad)+aHt./(sin(elRad)+bHt./(sin(elRad)+cHt)));

% adjust by 1000 because h should be in KM
dm = (1./sin(elRad)-mht).*h/1000;

Md = Md0+dm;

end















