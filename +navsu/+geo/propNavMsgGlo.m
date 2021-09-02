function pos = propNavMsgGlo(eph, prn, GPSweek, GPSsec)



WARNING_ENABLE = false;

% mu = 3.986005e14; % WGS-84 value in m^3/s^2
% OMEGA_DOTe = 7.2921151467e-5; % WGS-84 value in rad/s

tmArray = GPSweek * 604800 + GPSsec;
tmArrayLen = length(tmArray);
pos = struct('x', NaN(tmArrayLen, 1), 'y', NaN(tmArrayLen, 1), ...
             'z', NaN(tmArrayLen, 1), 'x_dot', NaN(tmArrayLen, 1), ...
             'y_dot', NaN(tmArrayLen, 1), 'z_dot', NaN(tmArrayLen, 1), ...
             'clock_bias', NaN(tmArrayLen, 1), 'clock_drift', NaN(tmArrayLen, 1), ...
             'clock_rel', NaN(tmArrayLen, 1), ...
             'accuracy', NaN(tmArrayLen, 1), 'health', NaN(tmArrayLen, 1), ...
             'IODC', NaN(tmArrayLen, 1), 't_m_toe', NaN(tmArrayLen, 1), ...
             'tslu', NaN(tmArrayLen, 1), 'toe_m_ttom', NaN(tmArrayLen, 1), ...
             'Fit_interval', NaN(tmArrayLen, 1), 'TGD', NaN(tmArrayLen,1), ...
             't_m_toc',NaN(tmArrayLen,1), 'freqNum',NaN(tmArrayLen,1));
% get all involved weekdays
weekDays = unique(eph.GPS_weekday(isfinite(eph.GPS_weekday)));

for loop = 1:tmArrayLen
    % find the most recent prior ephemeris
    idx = find(eph.PRN == prn(loop));
    if isempty(idx) || isnan(tmArray(loop))
        continue
    end
    
    % Finding most recent prior ephemeris
    t = tmArray(loop) ...
        - eph.GPS_week_num(idx) * 604800 ...
        - eph.GPS_weekday(idx)*86400*0 ...
        - eph.tk(idx);
    
    t(t < 0) = Inf;
    [~,tdx] = min(t);
    I = idx(tdx);
    
    % which ephemeris day is it (can have multiple)
    Iday = weekDays == eph.GPS_weekday(I);
    if sum(Iday) ~= 1 || length(Iday) > length(eph.leapSecond)
        % couldn't find a day or doesn't match number of leap seconds
        Iday = 1;
    end
    
    %% Compute position and velocity at desired time
    % Propagate ECEF position and velocity to time of interest
    posi = eph.P(I,:)'/1000;
    vel = eph.V(I,:)'/1000;
    acc = eph.A(I,:)'/1000;
    acc(abs(acc) > 1e-5) = acc(abs(acc) > 1e-5)*1e-10;
%     toe = mod(eph.ToE(I),7*86400);
    tfinal = tmArray(loop) - eph.ToE(I) - eph.leapSecond(Iday);
    
    
    ho = sign(tfinal)*60;
    if ho == 0
        ho = 1;
    end
    
    nsteps = floor(tfinal/ho);
    h = ones(nsteps,1)*ho;
    if mod(tfinal,ho) ~=0
        h = [h; rem(tfinal,ho)];
    end
    
    % Earth rotation rate (rad/s)
    constants.we = 0.7292115e-4;
    % Earth gravitational constant
    constants.mu = 398600.44;
    % Second zonal coefficient C20 (-J2)
    constants.C20 = -1082.63e-6;
    % Earth radius
    constants.ae  = 6378.136;
    
    state = [posi; vel];
    % Propagate steps forward/backward using Runge-Kutta 4 integration
    for i = 1:length(h)
        k1 = orbitDifEqGlonass(state          ,acc,constants);
        k2 = orbitDifEqGlonass(state+h(i)*k1/2,acc,constants);
        k3 = orbitDifEqGlonass(state+h(i)*k2/2,acc,constants);
        k4 = orbitDifEqGlonass(state+h(i)*k3  ,acc,constants);
        state = state + h(i)/6*(k1+2*k2+2*k3+k4);
    end
    posi = state(1:6);
    
    % Translate position back to ITRF2000 if necessary
    % GLONASS used reference frame PZ-90.02 before December 31, 2013, 12:00
    % UTC, which requires frame shift.
    if tmArray(loop) < 1072526400
        posi(1:3) = posi(1:3) + [-0.36 0.08 0.18]'/1000;
    end
    
    pos.x(loop) = posi(1)*1000;
    pos.y(loop) = posi(2)*1000;
    pos.z(loop) = posi(3)*1000;
    pos.x_dot(loop) = posi(4)*1000;
    pos.y_dot(loop) = posi(5)*1000;
    pos.z_dot(loop) = posi(6)*1000;
    
    %% Propagate clock
    pos.clock_bias(loop) =  eph.clock_bias(I)+tfinal*eph.clock_drift(I);
    pos.clock_drift(loop) = eph.clock_drift(I);
    %% Relativisitic effect
    %     pos.clock_rel(loop) = -4.442807633e-10*e*eph.sqrtA(I)*sin(E);
    %     pos.TGD(loop) = eph.TGD(I);
    pos.t_m_toc(loop) = tfinal;
    %     pos.accuracy(loop) = eph.accuracy(I);
    pos.health(loop) = eph.health(I);
    pos.IODC(loop) = eph.ToE(I);
    %     pos.AoD(loop) = tmArray(loop)- eph.GPS_week_num(I) * 604800 - eph.GPS_weekday(I)*86400- eph.ToE(I);
    pos.t_m_toe(loop) = tmArray(loop)- eph.ToE(I);
    
    pos.tslu(loop) = tmArray(loop) ...
                     - eph.GPS_week_num(I) * 604800 ...
                     - eph.GPS_weekday(I) * 86400 ...
                     - eph.tk(I) ...
                     + eph.AoOper(I) * 86400;
                 
    pos.tslu(loop) = eph.AoOper(I)*86400;
    
    pos.toe_m_ttom(loop) = eph.ToE(I) ...
                           - eph.GPS_week_num(I) * 604800 ...
                           - eph.GPS_weekday(I) * 86400*0 ...
                           - eph.tk(I);

    pos.Fit_interval(loop) = 900;
    pos.freqNum(loop) = eph.freqNum(I);
end


end


function dstate = orbitDifEqGlonass(state,acc,constants)

% Unpack input vectors
x  = state(1);
y  = state(2);
z  = state(3);
vx = state(4);
vy = state(5);
vz = state(6);
ax = acc(1);
ay = acc(2);
az = acc(3);

mu = constants.mu;
C20 = constants.C20;
re  = constants.ae;
we = constants.we;

r = sqrt(x^2+y^2+z^2);
p = re/r;
xb = x/r;
yb = y/r;
zb = z/r;
mub = mu/r^2;

dvx = -mub*xb+3/2*C20*mub*xb*p^2*(1-5*zb^2)+ax+we^2*x+2*we*vy;
dvy = -mub*yb+3/2*C20*mub*yb*p^2*(1-5*zb^2)+ay+we^2*y-2*we*vx;
dvz = -mub*zb+3/2*C20*mub*zb*p^2*(3-5*zb^2)+az;

dstate = [vx vy vz dvx dvy dvz]';

end























