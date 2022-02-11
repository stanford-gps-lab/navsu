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
tmN = length(tmArray);
pos = struct('x', NaN(tmN, 1), 'y', NaN(tmN, 1), 'z', NaN(tmN, 1), ...
    'x_dot', NaN(tmN, 1), 'y_dot', NaN(tmN, 1), 'z_dot', NaN(tmN, 1), ...
    'clock_bias', NaN(tmN, 1), 'clock_drift', NaN(tmN, 1),'clock_rel', NaN(tmN, 1), ...
    'accuracy', NaN(tmN, 1), 'health', NaN(tmN, 1), ...
	'IODC', NaN(tmN, 1), 't_m_toe', NaN(tmN, 1), ...
    'tslu', NaN(tmN, 1), 'toe_m_ttom', NaN(tmN, 1), ...
    'AoD',NaN(tmN,1),...
    'Fit_interval', NaN(tmN, 1),'TGD',NaN(tmN,1),'t_m_toc',NaN(tmN,1));


[iEph, iLastU] = deal(zeros(tmN, 1));

% prnLast = 0;
for loop = 1:tmN
    % find the most recent prior ephemeris for this PRN
    idx = find(eph.PRN == prn(loop));

    if isempty(idx)
        continue
    end
    
    if all(isfinite(eph.TTOM(idx)))
        % use latest message received
        t = tmArray(loop) - eph.GPS_week_num(idx) * 604800 - eph.TTOM(idx);
        t(t < 0) = Inf;
    else
        % use eph closest to current time
        t = tmArray(loop) - eph.GPS_week_num(idx) * 604800 - eph.Toe(idx);
    end
    if ~any(isfinite(t))
        % possibly we set everything to inf!
        continue
    end
    [~, tdx] = min(abs(t));
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
                            < 14400
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

% only work with prns I have an ephemeris for
haveEph = iEph > 0;
if ~any(haveEph)
    % don't have any ephemeris for any of the satellites!
    return
end
iEph = iEph(haveEph);

%% check for GAL I/NAV vs. F/NAV
if strcmp(constellation, 'GAL')
    %  Parsing data source
    binData = dec2base(eph.codes_on_L2(iEph), 2, 10);
    
    % very naive. does not check for duplicates in fields!
    galEphSource = zeros(size(binData, 1), 1);
    galEphSource(binData(:, end) == '1')   = 1; % I/NAV E1-B
    galEphSource(binData(:, end-1) == '1') = 2; % F/NAV E5a-I
    galEphSource(binData(:, end-2) == '1') = 3; % I/NAV E5b-I
    
    galClkRef = zeros(size(binData, 1), 1);
    galClkRef(binData(:, end-8) == '1')   = 1; % af0-af2, Toc are for E5a,E1
    galClkRef(binData(:, end-9) == '1')   = 2; % af0-af2, Toc are for E5b,E1

    % check for consistency
    if any((galEphSource == 2 & galClkRef == 2) ...
         | (galEphSource == 3 & galClkRef == 1))
        error('Inconsistent GAL ephemeris data source bits.');
    end

    % mark entries that are for E5b, E1 clock
    galE5b = galClkRef == 2;
end


%% compute SV position

% get time since ephemeris
tsE = tmArray(haveEph) ...
      - eph.GPS_week_num(iEph) * 604800 ...
      - eph.Toe(iEph);
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

% propagate longitude of ascending node
OMEGA = eph.OMEGA(iEph) ...
      + (eph.OMEGA_DOT(iEph) - OMEGA_DOTe) .* tsE ...
      - OMEGA_DOTe * eph.Toe(iEph);

% corrected longitude of ascending node for some BDS satellites
if strcmp(constellation, 'BDS')
    b2c = prn(haveEph) <= 5;
    OMEGA(b2c) = OMEGA(b2c) - eph.OMEGA_DOT(iEph(b2c)).*tsE(b2c);
end

% compute SV position in ECEF
pos.x(haveEph) = xo .* cos(OMEGA) - yo .* cos(uri(:, 3)) .* sin(OMEGA);
pos.y(haveEph) = xo .* sin(OMEGA) + yo .* cos(uri(:, 3)) .* cos(OMEGA);
pos.z(haveEph) = yo .* sin(uri(:, 3));

if strcmp(constellation, 'BDS') && any(b2c)
    % need to further rotate Beidou GEOs
    phiTemp = -5*pi/180;
    RxTemp = [1 0 0; 0 cos(phiTemp) sin(phiTemp); 0 -sin(phiTemp) cos(phiTemp)];
    haveEphInt = find(haveEph);
    for i2c = find(b2c)'
        phiTemp = OMEGA_DOTe*tsE(i2c);
        RzTemp = [cos(phiTemp) sin(phiTemp) 0; -sin(phiTemp) cos(phiTemp) 0; 0 0 1];
        
        i2cE = haveEphInt(i2c);
        posTemp = RzTemp*RxTemp*[pos.x(i2cE) pos.y(i2cE) pos.z(i2cE)]';
        pos.x(i2cE) = posTemp(1);
        pos.y(i2cE) = posTemp(2);
        pos.z(i2cE) = posTemp(3);
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
pos.x_dot(haveEph) = xo_dot .* cos(OMEGA) - xo .* sin(OMEGA) .* OMEGA_dot ...
                     - yo_dot .* cos(uri(:, 3)) .* sin(OMEGA) ...
                     + yo .* sin(uri(:, 3)) .* sin(OMEGA) .* uri_dot(:, 3) ...
                     - yo .* cos(uri(:, 3)) .* cos(OMEGA) .* OMEGA_dot;
pos.y_dot(haveEph) = xo_dot .* sin(OMEGA) + xo .* cos(OMEGA) .* OMEGA_dot ...
                     + yo_dot .* cos(uri(:, 3)) .* cos(OMEGA) ...
                     - yo .* sin(uri(:, 3)) .* cos(OMEGA) .* uri_dot(:, 3) ...
                     - yo .* cos(uri(:, 3)) .* sin(OMEGA) .* OMEGA_dot;                      
pos.z_dot(haveEph) = yo_dot .* sin(uri(:, 3)) + yo .* cos(uri(:, 3)) .* uri_dot(:, 3);    

%% compute SV clock bias and drift
toC = tmArray(haveEph) ...
      - eph.GPS_week_num(iEph) * 604800 ...
      - eph.Toc(iEph); % time from clock reference epoch
pos.clock_drift(haveEph) = eph.clock_drift(iEph) + eph.clock_drift_rate(iEph) .* toC;

pos.clock_bias(haveEph) = eph.clock_bias(iEph) + pos.clock_drift(haveEph) .* toC;

if dualFreq 
    % Beidou dual frequency TGD offsets
    kf = (1561.098/1207.14)^2;
    tgdOffset = (eph.TGD2(iEph)-kf*eph.TGD(iEph)) ./ (1-kf);
    pos.clock_bias(haveEph) = pos.clock_bias(haveEph) - tgdOffset;
end

%% Relativisitic effect, other elements
pos.clock_rel(haveEph) = -4.442807633e-10 * e .* eph.sqrtA(iEph) .* sin(E);
if dualFreq
    pos.TGD(haveEph) = tgdOffset;
else
    pos.TGD(haveEph) = eph.TGD(iEph);
end
% check for Galileo E5b/E1 entries. Then group delay is stored in TGD2
% entry. See GAL ICD 5.1.5 and navsu.readfiles.loadRinexNav
if strcmp(constellation, 'GAL') && any(galE5b)
    ephIdx = find(haveEph);
    pos.TGD(ephIdx(galE5b)) = eph.TGD2(iEph(galE5b));
end
    
pos.t_m_toc(haveEph)     = toC;
pos.accuracy(haveEph)    = eph.accuracy(iEph);
pos.health(haveEph)      = eph.health(iEph);
pos.IODC(haveEph)        = eph.IODC(iEph);
pos.t_m_toe(haveEph)     = tsE;
pos.AoD(haveEph)         = tmArray(haveEph) ...
                           - eph.GPS_week_num(iEph) * 604800 ...
                           - eph.Toe(iEph); % what is this??

switch constellation 
    case 'GPS'
        pos.tslu(haveEph) = tmArray(haveEph) ...
                            - eph.GPS_week_num(iLastU(haveEph)) * 604800 ...
                            - eph.TTOM(iLastU(haveEph));

    case 'BDS'
        pos.tslu(haveEph) = aodi;
end
pos.toe_m_ttom(haveEph)  = eph.Toe(iEph) - eph.TTOM(iEph);
pos.Fit_interval(haveEph) = eph.Fit_interval(iEph);
    
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
