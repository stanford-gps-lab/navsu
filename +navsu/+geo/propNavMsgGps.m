function pos = propNavMsgGps(eph, prn, GPSweek, GPSsec, constellation,dualFreq)
%Propagate navigation broadcast for constellations that propagate like GPS.
% Accepts inputs for one constellation at a time.
% 
% propNavMsgGps(eph, prn, GPSweek, GPSsec, constellation,dualFreq)
%   
% INPUTS:
%   eph     - broadcast ephemeris structure for the constellation
%   prn     - Nx1 vector of desired PRN to propagate
%   GPSweek - Nx1 vector of GPS week numbers of propagation times
%   GPSsec  - Nx1 vector of GPS second of week of propagation times
%   constellation   - string identifying constellation to propagate. Has to
%                   be 'GPS', 'GAL' or 'BDS'. Defaults to 'GPS'.
%   dualFreq        - optional extra input flag for TGD adjustment for BDS

% It's probably easier just to call propNavMsg- :)



% Default to GPS
if nargin < 5
    constellation = 'GPS';
end

if nargin < 6
    % Flag to do dual frequency TGD adjustment for Beidou (default to
    % doing dual frequency adjustment for BDS)
    if strcmp(constellation,'BDS')
        dualFreq = 1;
    else
        dualFreq = 0;
    end
    
end

WARNING_ENABLE = false;

mu = 3.986005e14; % WGS-84 value in m^3/s^2
OMEGA_DOTe = 7.2921151467e-5; % WGS-84 value in rad/s

if strcmp(constellation,'BDS')
    OMEGA_DOTe = 7.292115e-5; 
end

tmArray = GPSweek * 604800 + GPSsec;
tmArrayLen = length(tmArray);
pos = struct('x', NaN(tmArrayLen, 1), 'y', NaN(tmArrayLen, 1), 'z', NaN(tmArrayLen, 1), ...
    'x_dot', NaN(tmArrayLen, 1), 'y_dot', NaN(tmArrayLen, 1), 'z_dot', NaN(tmArrayLen, 1), ...
    'clock_bias', NaN(tmArrayLen, 1), 'clock_drift', NaN(tmArrayLen, 1),'clock_rel', NaN(tmArrayLen, 1), ...
    'accuracy', NaN(tmArrayLen, 1), 'health', NaN(tmArrayLen, 1), ...
	'IODC', NaN(tmArrayLen, 1), 't_m_toe', NaN(tmArrayLen, 1), ...
    'tslu', NaN(tmArrayLen, 1), 'toe_m_ttom', NaN(tmArrayLen, 1), ...
    'AoD',NaN(tmArrayLen,1),...
    'Fit_interval', NaN(tmArrayLen, 1),'TGD',NaN(tmArrayLen,1),'t_m_toc',NaN(tmArrayLen,1));

if strcmp(constellation,'GAL')
   %  Parsing data source
   binData = repmat(' ',length(eph.codes_on_L2),22);
   binData(:,2:2:end) = dec2bin(eph.codes_on_L2,11);
   binData = str2num(binData);
   
   % very naive. does not check for duplicates in fields!
   galEphSource = zeros(length(eph.codes_on_L2),1);
   galEphSource(find(binData(:,end)))   = 1; % I/NAV E1-B
   galEphSource(find(binData(:,end-1))) = 2; % F/NAV E5a-I
   galEphSource(find(binData(:,end-2))) = 3; % I/NAV E5b-I
   
   galClkRef = zeros(length(eph.codes_on_L2),1);
   galClkRef(find(binData(:,end-8)))   = 1; % af0-af2, Toc are for E5a,E1
   galClkRef(find(binData(:,end-9)))   = 2; % af0-af2, Toc are for E5b,E1

end

[iEph, iLastU] = deal(zeros(tmArrayLen, 1));

% prnLast = 0;
for loop = 1:tmArrayLen
    % find the most recent prior ephemeris
%     if prnLast ~= prn(loop)
        switch constellation
            case 'GPS'
                idx = find(eph.PRN == prn(loop));
            case 'GAL'
                % Need to use af0 and af1 for E5a (CODE using E5a, not b)
                idx = find(eph.PRN == prn(loop));% & galClkRef == 1 & galEphSource == 2);
                
            case 'BDS'
                idx = find(eph.PRN == prn(loop));
        end
%     end
    if isempty(idx)
        continue
    end
    
    t = tmArray(loop) - eph.GPS_week_num(idx) * 604800 - eph.TTOM(idx);
    t(t < 0) = Inf;
    [~,tdx] = min(t);
    I = idx(tdx);
    
    % look for new uploads (or check secondary health set)
    switch constellation 
        case 'GPS'
            % find the most recent ephemeris associated with an upload
            jdx = find(eph.PRN == prn(loop) & mod(eph.Toe,7200) ~= 0);
            if ~isempty(jdx)
                % after an upload two ephemerides may have the offset toe. Only use
                % the first one
                k = 2;
                while k <= length(jdx)
                    if mod(eph.Toe(jdx(k)),7200) == mod(eph.Toe(jdx(k-1)),7200) ...
                            && abs(eph.TTOM(jdx(k)) + eph.GPS_week_num(jdx(k)) - ...
                            eph.TTOM(jdx(k-1)) - eph.GPS_week_num(jdx(k-1))) ...
                            < 14400;
                        if eph.TTOM(jdx(k)) + eph.GPS_week_num(jdx(k)) >= ...
                                eph.TTOM(jdx(k-1)) - eph.GPS_week_num(jdx(k-1))
                            jdx(k) = [];
                        else
                            jdx(k-1) = [];
                        end
                    end
                    k = k+1;
                end
                tmp = tmArray(loop) - eph.GPS_week_num(jdx) * 604800 - eph.TTOM(jdx);
                tmp(tmp < 0) = max(tmp) - tmp(tmp < 0);
                %     t(t < 60) = max(t) - t(t < 60);
                [~, K] = min(tmp(end:-1:1));
                K = jdx(length(tmp) + 1 - K);
            else
                [~, tmp] = min(eph.GPS_week_num(idx) * 604800 + eph.TTOM(idx));
                K = idx(tmp);
            end            
            
        case 'GAL'
            % no known data to suggest uploads 
%             t = tmArray(loop) - eph.GPS_week_num(idx) * 604800 - eph.TTOM(idx2);
%             t(t < 0) = Inf;
%             [~,tdx] = min(t);
%             I = idx(tdx2);
            
%             if prn(loop) == 30 && tmArray(
                
            
        case 'BDS'
            % pull age of data straight from nav data (need to translate)
%             aodi = eph.aodClock(I)*3600;
            aodi = 0;
            
    end
    % check health
    if WARNING_ENABLE && eph.health(I)
        fprintf(2, 'Warning: SV health %g for the ephemeris No.%d\n', eph.health(I), I);
    end
    
    % save the important indices for this satellite
    iEph(loop) = I;
    if strcmp(constellation, 'GPS')
        iLastU(loop) = K;
    end
    
end

% the remaining calculations can be done vectorized!

% time since ephemeris for each SV
tsE = tmArray - eph.GPS_week_num(iEph) * 604800 - eph.Toe(iEph);
if strcmp(constellation,'BDS')
    tsE = tsE - 14; % need to adjust for Beidou timescale leapseconds
end

% check fit interval
if WARNING_ENABLE 
    for iWarn = find(tsE > eph.Fit_interval(iEph) * 3600)'
        fprintf(2, 'Warning: t = %g (second) exceeds the Fit Interval %g (hour)\n', ...
            tsE(iWarn), eph.Fit_interval(iWarn));
    end
end

%% compute SV position

A = eph.sqrtA(iEph).^2; % semimajor axis  
n0 = sqrt(mu ./ A.^3); % mean motion in rad/s

if isfield(eph,'aDot')
    A = A + eph.aDot(iEph) .* tsE; % correct A with Adot if CNAV
end

n = n0 + eph.Delta_n(iEph); % corrected mean motion
M = eph.M0(iEph) + n .* tsE; % mean anomaly
e = eph.e(iEph); % eccentricity
E = KeplerEq(M, e); % solve Kepler's eq for eccentric anomaly
nu = atan2(sqrt(1 - e.^2) .* sin(E), cos(E) - e); % true anomaly
PHI = nu + eph.omega(iEph); % argument of latitude

% correct argument of latitude, radius, and inclination      
d_uri = [dot([eph.Cus(iEph) eph.Cuc(iEph)], [sin(2*PHI) cos(2*PHI)], 2) ...
         dot([eph.Crs(iEph) eph.Crc(iEph)], [sin(2*PHI) cos(2*PHI)], 2) ...
         dot([eph.Cis(iEph) eph.Cic(iEph)], [sin(2*PHI) cos(2*PHI)], 2)];

uri = d_uri + [PHI, A .* (1 - e .* cos(E)), eph.i0(iEph) + eph.I_DOT(iEph).*tsE];

xo = uri(:, 2) .* cos(uri(:, 1)); % SV x-position in orbital plane
yo = uri(:, 2) .* sin(uri(:, 1)); % SV y-position in orbital plane

% corrected longitude of ascending node for some BDS satellites

OMEGA = eph.OMEGA(iEph) ...
      + (eph.OMEGA_DOT(iEph) - OMEGA_DOTe) .* tsE...
      - OMEGA_DOTe * eph.Toe(iEph);

if strcmp(constellation, 'BDS')
    b2c = prn <= 5;
    OMEGA(b2c) = OMEGA(b2c) - eph.OMEGA_DOT(iEph(b2c)).*tsE(b2c);
end

% compute SV position in ECEF
pos.x = xo .* cos(OMEGA) - yo .* cos(uri(:, 3)) .* sin(OMEGA);
pos.y = xo .* sin(OMEGA) + yo .* cos(uri(:, 3)) .* cos(OMEGA);
pos.z = yo .* sin(uri(:, 3));

if strcmp(constellation, 'BDS') && any(b2c)
    % need to further rotate Beidou GEOs
    phiTemp = -5*pi/180;
    RxTemp = [1 0 0; 0 cos(phiTemp) sin(phiTemp); 0 -sin(phiTemp) cos(phiTemp)];
    for i2c = find(b2c)'
        phiTemp = OMEGA_DOTe*tsE(i2c);
        RzTemp = [cos(phiTemp) sin(phiTemp) 0; -sin(phiTemp) cos(phiTemp) 0; 0 0 1];

        posTemp = RzTemp*RxTemp*[pos.x(i2c) pos.y(i2c) pos.z(i2c)]';
        pos.x(i2c) = posTemp(1);
        pos.y(i2c) = posTemp(2);
        pos.z(i2c) = posTemp(3);
    end
end

%% compute SV velocity
M_dot = n; % rate of mean anomaly
E_dot = M_dot ./ (1 - e .* cos(E)); % rate of eccentric anomaly
nu_dot = E_dot .* sqrt(1 - e.^2) ./ (1 - e .* cos(E)); % rate of true anomaly
PHI_dot = nu_dot; % rate of argument of latitude

% rates of argument of latitude, radius, and inclination             
d_uri_dot = [dot([eph.Cus(iEph) eph.Cuc(iEph)], [cos(2*PHI), -sin(2*PHI)], 2) ...
             dot([eph.Crs(iEph) eph.Crc(iEph)], [cos(2*PHI), -sin(2*PHI)], 2) ...
             dot([eph.Cis(iEph) eph.Cic(iEph)], [cos(2*PHI), -sin(2*PHI)], 2)] ...
           .* (2*PHI_dot);

uri_dot = d_uri_dot + [PHI_dot, A .* e .* sin(E) .* E_dot, eph.I_DOT(iEph)];

% SV x-velocity in orbital plane
xo_dot = uri_dot(:, 2) .* cos(uri(:, 1)) - uri(:, 2) .* sin(uri(:, 1)) .* uri_dot(:, 1); 
% SV y-velocity in orbital plane
yo_dot = uri_dot(:, 2) .* sin(uri(:, 1)) + uri(:, 2) .* cos(uri(:, 1)) .* uri_dot(:, 1); 
% rate of longitude of ascending node
OMEGA_dot = eph.OMEGA_DOT(iEph) - OMEGA_DOTe;   

% compute SV velocity in ECEF
pos.x_dot = xo_dot .* cos(OMEGA) - xo .* sin(OMEGA) .* OMEGA_dot ...
          - yo_dot .* cos(uri(:, 3)) .* sin(OMEGA) ...
          + yo .* sin(uri(:, 3)) .* sin(OMEGA) .* uri_dot(:, 3) ...
          - yo .* cos(uri(:, 3)) .* cos(OMEGA) .* OMEGA_dot;
pos.y_dot = xo_dot .* sin(OMEGA) + xo .* cos(OMEGA) .* OMEGA_dot ...
          + yo_dot .* cos(uri(:, 3)) .* cos(OMEGA) ...
          - yo .* sin(uri(:, 3)) .* cos(OMEGA) .* uri_dot(:, 3) ...
          - yo .* cos(uri(:, 3)) .* sin(OMEGA) .* OMEGA_dot;                      
pos.z_dot = yo_dot .* sin(uri(:, 3)) + yo .* cos(uri(:, 3)) .* uri_dot(:, 3);    

%% compute SV clock bias and drift
toC = tmArray - eph.GPS_week_num(iEph) * 604800 - eph.Toc(iEph); % time from clock reference epoch
pos.clock_drift = eph.clock_drift(iEph) + eph.clock_drift_rate(iEph) .* toC;

pos.clock_bias = eph.clock_bias(iEph) + pos.clock_drift .* toC;

if dualFreq 
    % Beidou dual frequency TGD offsets
    kf = (1561.098/1207.14)^2;
    tgdOffset = (eph.TGD2(iEph)-kf*eph.TGD(iEph)) ./ (1-kf);
    pos.clock_bias = pos.clock_bias - tgdOffset;
end

%% Relativisitic effect, other elements
pos.clock_rel = -4.442807633e-10 * e .* eph.sqrtA(iEph) .* sin(E);
if dualFreq
    pos.TGD = tgdOffset;
else
    pos.TGD = eph.TGD(iEph);
end
pos.t_m_toc     = toC;
pos.accuracy    = eph.accuracy(iEph);
pos.health      = eph.health(iEph);
pos.IODC        = eph.IODC(iEph);
pos.t_m_toe     = tsE;
pos.AoD         = tmArray - eph.GPS_week_num(iEph) * 604800 - eph.Toe(iEph); % what is this??

switch constellation 
    case 'GPS'
        pos.tslu = tmArray - eph.GPS_week_num(iLastU) * 604800 - eph.TTOM(iLastU);

    case 'BDS'
        pos.tslu = aodi;
end
pos.toe_m_ttom  = eph.Toe(iEph) - eph.TTOM(iEph);
pos.Fit_interval = eph.Fit_interval(iEph);
    
end

% end

function E = KeplerEq(M0, e)
% E = KeplerEq(M, e);
%
% Function (vectorized) to iteratively solve Kepler's 
% equation for eccentric anomaly.
%
% Input:
%   M - mean anomaly (rad) (vector)
%   e - eccentricity (dimensionless) (nxm) or (1x1)
% Output:
%   E - eccentric anomaly (rad) (nxm)
%

% wrapToPi = @(x) mod(x + pi, 2 * pi) - pi;

e = abs(e);
M =  mod(M0 + pi, 2 * pi) - pi;
Ms = sign(M);
Mr = M0 - M;
M = abs(M);

% set convergence criterion
tol = 1e-15;      % convergence criterion
max_iter = 10;    % max iteration 

m = 1;            % iteration counter initial value
dE = ones(size(e));         % set initial value of delta-E to fail convergence
nConv = abs(dE) > tol; % index of satellites that have not yet converged

E = M + (e .* sin(M)) ./ (1 - sin(M + e) + sin(M));    % set initial value of E 

while any(nConv) && (m <= max_iter)
  dE = (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
  E(nConv) = E(nConv) - dE(nConv);         % update E
  m = m + 1;          % increment iteration numbers
  nConv = abs(dE) > tol; % check again for satellites that have not converged

  if m > max_iter     % check against maximum iterations allowed
    fprintf(2, 'Warning message from KEPLR_EQ ...\n')
    fprintf(2, 'Maximum iterations exceeded in solution of Kepler''s Equation.\n')
    fprintf(2, 'Results may be invalid.\n\n')
    return
  end % if

end % while

E = Mr + E .* Ms;

end
