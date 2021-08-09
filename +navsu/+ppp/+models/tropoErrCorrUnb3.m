function [tropo_corr, Mw,tropDataExtra] = tropoErrCorrUnb3(el, h, lat,doy)

% tropo_corr = zeros(size(el));

% 
latTable = [15 30 45 60 75];
% interpolation tables
       %   p_o    T_o    e_o   B_o   lam_o
table1 = [1013.25 299.65 26.31 6.30e-3 2.77;
          1017.25 294.15 21.79 6.05e-3 3.15;
          1015.75 283.15 11.66 5.58e-3 2.57;
          1011.75 272.15 6.78  5.39e-3 1.81;
          1013.00 263.65 4.11  4.53e-3 1.55];
       
      % delta of above values in table 1
table2 = [ 0.00  0.00 0.00 0.00e-3 0.00
          -3.75  7.00 8.85 0.25e-3 0.33
          -2.25 11.00 7.24 0.32e-3 0.46
          -1.75 15.00 5.36 0.81e-3 0.74
          -0.50 14.50 3.39 0.62e-3 0.30   ];
      
% Interpolating all values 
% get lat values for interpolation
absLat = max(min(abs(lat), latTable(end)), latTable(1));

metValue0 = interp1(latTable', table1, absLat);
dmetValue = interp1(latTable', table2, absLat);

Dmin = 28*ones(size(el));
Dmin(lat < 0) = 211;
metValue = metValue0-dmetValue.*cos((2*pi*(doy-Dmin))/365.25);

B = metValue(:,4);
T = metValue(:,2);
p = metValue(:,1);
e = metValue(:,3);
lam = metValue(:,5);

% Constants
k1 = 77.604;  % K/mbar
k2 = 382000;  % K^2/mbar
Rd = 287.054; % J/kg/K
gm = 9.784;   % m/s^2
g  = 9.80665; %m/s^2

d_dry = (1-B.*h./T).^(g./(Rd.*B)).*((10^-6*k1*Rd*p)/gm);
% d_wet = (1-B.*h./T).^((lam+1).*g./(Rd.*B)).*((10^-6*k1*Rd*p.*e)./(T.*(gm.*(lam+1)-B.*Rd)));
d_wet = (1-B.*h./T).^((lam+1).*g./(Rd.*B)-1).*((10^-6*k2*Rd.*e)./(T.*(gm.*(lam+1)-B.*Rd)));

% m = 1.001./sqrt(0.002001 + sin(el*pi/180).^2);

[Mw, Md] = navsu.ppp.models.mappingOfNiell(el, h, lat,doy);

% Mw = m(:,1);
% Md = m(:,2);

tropo_corr = Md.*d_dry+Mw.*d_wet;
% tropo_corr = tropo_corr;

% Mw = Mw;

tropDataExtra.trototSave = d_dry;
tropDataExtra.gmfwSave   = Mw;
tropDataExtra.tzd        = d_dry(1)+d_wet(1);
tropDataExtra.ddry       = d_dry(1);
tropDataExtra.dwet       = d_wet(1);

% remove bad data...
tropo_corr(tropo_corr > 1e3) = 0;
if ~isreal(tropo_corr)
    tropo_corr = zeros(size(tropo_corr));
end

end
