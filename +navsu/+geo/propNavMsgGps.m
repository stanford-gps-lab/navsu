function pos = propNavMsgGps(eph, prn, GPSweek, GPSsec, constellation,dualFreq)

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

prnLast = 0;
for loop = 1:tmArrayLen
    % find the most recent prior ephemeris
%     if prnLast ~= prn(loop)
        switch constellation
            case 'GPS'
                idx = find(eph.PRN == prn(loop));
            case 'GAL'
                % Need to use af0 and af1 for E5a (CODE using E5a, not b)
                idx = find(eph.PRN == prn(loop) & galClkRef == 1 & galEphSource == 2);
                
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
    
    A = eph.sqrtA(I)^2; % semimajor axis  
    n0 = sqrt(mu / A^3); % mean motion in rad/s
    if strcmp(constellation,'BDS')
        % need to adjust for Beidou timescale leapseconds
        t = tmArray(loop) - eph.GPS_week_num(I) * 604800 - eph.Toe(I)-14; % time from ephemeris reference epoch
    else
        t = tmArray(loop) - eph.GPS_week_num(I) * 604800 - eph.Toe(I); % time from ephemeris reference epoch
    end
%     t = mod(t + 302400, 604800) - 302400; % !!
    if isfield(eph,'aDot')
        A = A + eph.aDot(I)*t; % correct A with Adot if CNAV
    end


    % check fit interval
    if WARNING_ENABLE && t > eph.Fit_interval(I) * 3600
        fprintf(2, 'Warning: t = %g (second) exceeds the Fit Interval %g (hour)\n', t, eph.Fit_interval(I));
    end
    
    n = n0 + eph.Delta_n(I); % corrected mean motion
    M = eph.M0(I) + n * t; % mean anomaly
    e = eph.e(I); % eccentricity
    E = KeplerEq(M, e); % solve Kepler's eq for eccentric anomaly
    nu = atan2(sqrt(1 - e.^2) * sin(E), cos(E) - e); % true anomaly
    PHI = nu + eph.omega(I); % argument of latitude
    
    % correct argument of latitude, radius, and inclination
    d_uri = [eph.Cus(I) eph.Cuc(I); eph.Crs(I) eph.Crc(I); eph.Cis(I) eph.Cic(I)] ...
        * [sin(2*PHI); cos(2*PHI)];
    uri = d_uri + [PHI; A .* (1 - e .* cos(E)); eph.i0(I) + eph.I_DOT(I)*t];
    
    xo = uri(2) .* cos(uri(1)); % SV x-position in orbital plane
    yo = uri(2) .* sin(uri(1)); % SV y-position in orbital plane
    
    % corrected longitude of ascending node
     if strcmp(constellation,'BDS') && prn(loop) <= 5
        OMEGA = eph.OMEGA(I) + (eph.OMEGA_DOT(I)) * t - OMEGA_DOTe * eph.Toe(I);
        
        % compute SV position in ECEF
        pos.x(loop) = xo * cos(OMEGA) - yo * cos(uri(3)) * sin(OMEGA);
        pos.y(loop) = xo * sin(OMEGA) + yo * cos(uri(3)) * cos(OMEGA);
        pos.z(loop) = yo * sin(uri(3));
        % need to further rotate Beidou GEOs
        phiTemp = -5*pi/180;
        RxTemp = [1 0 0; 0 cos(phiTemp) sin(phiTemp); 0 -sin(phiTemp) cos(phiTemp)];
        phiTemp = OMEGA_DOTe*t;
        RzTemp = [cos(phiTemp) sin(phiTemp) 0; -sin(phiTemp) cos(phiTemp) 0; 0 0 1];
        
        posTemp = RzTemp*RxTemp*[pos.x(loop) pos.y(loop) pos.z(loop)]';
        pos.x(loop) = posTemp(1);
        pos.y(loop) = posTemp(2);
        pos.z(loop) = posTemp(3);
    else
        OMEGA = eph.OMEGA(I) + (eph.OMEGA_DOT(I) - OMEGA_DOTe) * t - OMEGA_DOTe * eph.Toe(I);
        
        % compute SV position in ECEF
        pos.x(loop) = xo * cos(OMEGA) - yo * cos(uri(3)) * sin(OMEGA);
        pos.y(loop) = xo * sin(OMEGA) + yo * cos(uri(3)) * cos(OMEGA);
        pos.z(loop) = yo * sin(uri(3));
    end

    
    M_dot = n; % rate of mean anomaly
    E_dot = M_dot ./ (1 - e .* cos(E)); % rate of eccentric anomaly
    nu_dot = E_dot .* sqrt(1 - e.^2) ./ (1 - e .* cos(E)); % rate of true anomaly
    PHI_dot = nu_dot; % rate of argument of latitude

    % rates of argument of latitude, radius, and inclination
    d_uri_dot = ([eph.Cus(I) eph.Cuc(I); eph.Crs(I) eph.Crc(I); ...
                   eph.Cis(I) eph.Cic(I)] ...
                 * [cos(2*PHI); -sin(2*PHI)])*(2*PHI_dot);
    uri_dot = d_uri_dot + [PHI_dot; A .* e .* sin(E) .* E_dot; eph.I_DOT(I)];
           
    % SV x-velocity in orbital plane
    xo_dot = uri_dot(2) .* cos(uri(1)) - uri(2) .* sin(uri(1)) .* uri_dot(1); 
    % SV y-velocity in orbital plane
    yo_dot = uri_dot(2) .* sin(uri(1)) + uri(2) .* cos(uri(1)) .* uri_dot(1); 
    % rate of longitude of ascending node
    OMEGA_dot = eph.OMEGA_DOT(I) - OMEGA_DOTe;   
    
    % compute SV velocity in ECEF
    pos.x_dot(loop) = xo_dot * cos(OMEGA) - xo * sin(OMEGA) * OMEGA_dot ...
                      - yo_dot * cos(uri(3)) * sin(OMEGA) ...
                      + yo * sin(uri(3)) * sin(OMEGA) * uri_dot(3) ...
                      - yo * cos(uri(3)) * cos(OMEGA) * OMEGA_dot;
    pos.y_dot(loop) = xo_dot * sin(OMEGA) + xo * cos(OMEGA) * OMEGA_dot ...
                      + yo_dot * cos(uri(3)) * cos(OMEGA) ...
                      - yo * sin(uri(3)) * cos(OMEGA) * uri_dot(3) ...
                      - yo * cos(uri(3)) * sin(OMEGA) * OMEGA_dot;                      
    pos.z_dot(loop) = yo_dot * sin(uri(3)) + yo * cos(uri(3)) * uri_dot(3);    

    % compute SV clock bias and drift
    t = tmArray(loop) - eph.GPS_week_num(I) * 604800 - eph.Toc(I); % time from clock reference epoch
    pos.clock_drift(loop) = eph.clock_drift(I) + eph.clock_drift_rate(I) * t;
    
    if dualFreq 
        % Beidou dual frequency TGD offsets
        kf = (1561.098/1207.14)^2;
        tgdOffset = (eph.TGD2(I)-kf*eph.TGD(I))./(1-kf);
        pos.clock_bias(loop) = eph.clock_bias(I) + pos.clock_drift(loop) * t-tgdOffset;
    else
        pos.clock_bias(loop) = eph.clock_bias(I) + pos.clock_drift(loop) * t;
    end
    % pos.clock_bias(loop) = eph.clock_bias(I) + pos.clock_drift(loop) * t - 4.442807633e-10 * e * eph.sqrtA(I) * sin(E);
    
    % Relativisitic effect
    pos.clock_rel(loop) = -4.442807633e-10*e*eph.sqrtA(I)*sin(E);
    if dualFreq
        pos.TGD(loop) = tgdOffset;
    else
        pos.TGD(loop) = eph.TGD(I);
    end
    pos.t_m_toc(loop) = t;
    pos.accuracy(loop) = eph.accuracy(I);
    pos.health(loop) = eph.health(I);
    pos.IODC(loop) = eph.IODC(I);
    pos.t_m_toe(loop) = tmArray(loop) - eph.GPS_week_num(I) * 604800 - eph.Toe(I);
    pos.AoD(loop) = tmArray(loop) - eph.GPS_week_num(I) * 604800 - eph.Toe(I);
    
    switch constellation 
        case 'GPS'
            pos.tslu(loop) = tmArray(loop) - eph.GPS_week_num(K) * 604800 - eph.TTOM(K);
            
        case 'BDS'
            pos.tslu(loop) = aodi;
    end
    pos.toe_m_ttom(loop) = eph.Toe(I) - eph.TTOM(I);
    pos.Fit_interval(loop) = eph.Fit_interval(I);
    
end

end

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
dE = 1.0;         % set initial value of delta-E to fail convergenge

E = M + (e .* sin(M)) ./ (1 - sin(M + e) + sin(M));    % set initial value of E 

while any(abs(dE) > tol) && (m <= max_iter)
  dE = (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
  E = E - dE;         % update E
  m = m + 1;          % increment iteration numbers

  if m > max_iter     % check against maximum iterations allowed
    fprintf(2, 'Warning message from KEPLR_EQ ...\n')
    fprintf(2, 'Maximum iterations exceeded in solution of Kepler''s Equation.\n')
    fprintf(2, 'Results may be invalid.\n\n')
    return
  end % if

end % while

E = Mr + E .* Ms;

end
