function timeUpdate(filter,epoch)
% PPP filter time update


% pull the position, velocity, and covariance off
pos   = filter.pos;
vel   = filter.vel;
cov  = filter.cov;


PARAMS = filter.PARAMS;

% Time since last GNSS kalman filter update
dtKf = epoch-filter.epochLastGnssUpdate;

filter.epochLastGnssUpdate = epoch;

% Constants (sone of these could be changed to inputs at a later date)

% inertial updates have not occurred - this is a normal KF time and
% measurement update
Phi_matrix = eye(size(cov));
% velocity -> position
Phi_matrix(filter.INDS_STATE.POS,filter.INDS_STATE.VEL) = eye(3) * dtKf;
% clock drift -> clock bias
Phi_matrix(filter.INDS_STATE.CLOCK_BIAS,filter.INDS_STATE.CLOCK_DRIFT) = dtKf;

% 2. Determine approximate system noise covariance matrix using (14.82)
Q              = zeros(size(cov));
Q(filter.INDS_STATE.ATTITUDE,filter.INDS_STATE.ATTITUDE) = ...
    eye(3) * PARAMS.Q.gyro_noise_PSD.^2 * dtKf;
Q(filter.INDS_STATE.POS,filter.INDS_STATE.POS)     =  ...
    eye(3) * PARAMS.Q.POS.^2 * dtKf;
Q(filter.INDS_STATE.VEL,filter.INDS_STATE.VEL)     =  ...
    eye(3) * PARAMS.Q.VEL.^2 * dtKf;
Q(filter.INDS_STATE.ACC_BIAS,filter.INDS_STATE.ACC_BIAS) = ...
    eye(3) * PARAMS.Q.ACC_BIAS.^2 * dtKf;
Q(filter.INDS_STATE.W_BIAS,filter.INDS_STATE.W_BIAS) = ...
    eye(3) * PARAMS.Q.W_BIAS.^2 * dtKf;
Q(filter.INDS_STATE.CLOCK_BIAS,filter.INDS_STATE.CLOCK_BIAS) = ...
    eye(length(filter.INDS_STATE.CLOCK_BIAS)).*(PARAMS.Q.RXB).^2 * dtKf;
Q(filter.INDS_STATE.CLOCK_DRIFT,filter.INDS_STATE.CLOCK_DRIFT) = ...
    eye(length(filter.INDS_STATE.CLOCK_BIAS)).*(PARAMS.Q.RXB).^2 * dtKf;

if PARAMS.states.RX_DCB
    Q(filter.INDS_STATE.RX_DCB.INDS,filter.INDS_STATE.RX_DCB.INDS)     =  ...
        eye(length(filter.INDS_STATE.RX_DCB.INDS)) * PARAMS.Q.RX_DCB.^2 * dtKf;
end

if PARAMS.states.trop
    Q(filter.INDS_STATE.TROP,filter.INDS_STATE.TROP)     =  ...
        eye(length(filter.INDS_STATE.TROP)) * PARAMS.Q.TROP.^2 * dtKf;
end

if PARAMS.states.iono && strcmp(PARAMS.states.ionoMode,'L1DELAYSTATE')
    indsTec = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
    Q(indsTec,indsTec) = eye(length(indsTec))*PARAMS.Q.L1_IONO.^2*dtKf;
end

if PARAMS.states.iono && strcmp(PARAMS.states.ionoMode,'TECSTATE')
    indsTec = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
    Q(indsTec,indsTec) = eye(length(indsTec))*PARAMS.Q.L1_IONO.^2*dtKf;
end

if PARAMS.states.RX_DCB_GLO
    inds_RX_DCB_GLO = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 3 & filter.INDS_STATE.FLEX_STATES_INFO(:,2) == 2);
    Q(inds_RX_DCB_GLO,inds_RX_DCB_GLO) = eye(length(inds_RX_DCB_GLO))*PARAMS.Q.RX_DCB_GLO.^2*dtKf;
end

if PARAMS.states.RX_DCB_GPS
    inds_RX_DCB_GPS = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 3 & ...
        filter.INDS_STATE.FLEX_STATES_INFO(:,2) == 1);
    Q(inds_RX_DCB_GPS,inds_RX_DCB_GPS) = eye(length(inds_RX_DCB_GPS))*PARAMS.Q.RX_DCB_GPS.^2*dtKf;
end

if PARAMS.states.MP_CODE
    inds_MP_CODE = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 4);
    
    Phi_matrix(inds_MP_CODE,inds_MP_CODE) = eye(length(inds_MP_CODE))*exp(-dtKf/PARAMS.other.TAU_MP_CODE);
    
    Q(inds_MP_CODE,inds_MP_CODE) = eye(length(inds_MP_CODE))*PARAMS.Q.MP_CODE.^2*dtKf;
end

if PARAMS.states.MP_CARR
    inds_MP_CARR = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 5);
    Phi_matrix(inds_MP_CARR,inds_MP_CARR) = eye(length(inds_MP_CARR))*exp(-dtKf/PARAMS.other.TAU_MP_CODE);
    Q(inds_MP_CARR,inds_MP_CARR) = eye(length(inds_MP_CARR))*PARAMS.Q.MP_CARR.^2*dtKf;
end

indsAmbs = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 1);
Q(indsAmbs,indsAmbs) = eye(length(indsAmbs))*PARAMS.Q.AMB.^2*dtKf;

if PARAMS.states.EPH
    indsEph = filter.INDS_STATE.FLEX_STATES(filter.INDS_STATE.FLEX_STATES_INFO(:,3) == 6);
    Q(indsEph,indsEph) = eye(length(indsEph))*PARAMS.Q.EPH.^2*dtKf;
end

% Propagate position
pos = pos+dtKf*vel;

% 3. Propagate state estimates using (3.14) noting that only the clock
% states are non-zero due to closed-loop correction.
x_est_propagated = zeros(size(filter.state));
x_est_propagated(filter.INDS_STATE.CLOCK_BIAS,1)   = filter.clockBias+dtKf*filter.clockDrift;
x_est_propagated(filter.INDS_STATE.CLOCK_DRIFT,1)  = filter.clockDrift;
x_est_propagated(filter.INDS_STATE.TROP,1)         = filter.state(filter.INDS_STATE.TROP);
x_est_propagated(filter.INDS_STATE.FLEX_STATES,1)  = filter.state(filter.INDS_STATE.FLEX_STATES);
x_est_propagated(filter.INDS_STATE.RX_DCB.INDS,1)  = filter.state(filter.INDS_STATE.RX_DCB.INDS);

% 4. Propagate state estimation error covariance matrix using (3.46)
cov_propagated = Phi_matrix * (cov + 0.5 * Q) *...
    Phi_matrix' + 0.5 * Q;


%% Put everything back in the filter
filter.pos = pos;
filter.vel = vel;
filter.clockBias = x_est_propagated(filter.INDS_STATE.CLOCK_BIAS,1);

filter.state = x_est_propagated;
filter.cov   = cov_propagated;





end