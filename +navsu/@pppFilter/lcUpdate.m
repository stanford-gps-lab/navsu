function lcUpdate(obj,epoch,posMeasi,velMeasi,PARAMS)

R_b_e = obj.R_b_e;
pos   = obj.pos;
vel   = obj.vel;
imuBiasStates = obj.imuBiasStates;
cov  = obj.cov;
nState = size(cov,1);

if isempty(obj.lastAccMeas)
    accMeasi = zeros(3,1);
else
    accMeasi = obj.lastAccMeas;
end
if isempty(obj.lastGyroMeas)
    gyroMeasi = zeros(3,1);
else
    gyroMeasi = obj.lastGyroMeas;
end

llhi = xyz2llh(pos');
latRad = llhi(1)*pi/180;

dtKf = epoch-obj.epochLastGnssUpdate;
obj.epochLastGnssUpdate = epoch;      
        
% Constants (sone of these could be changed to inputs at a later date)
c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity                
  
% Skew symmetric matrix of Earth rate
Omega_ie = crossProdMatrix([0,0,omega_ie]);          
            
             
% 1. Determine transition matrix using (14.50) (first-order approx)
Phi_matrix = eye(nState);
Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) - Omega_ie * dtKf;
Phi_matrix(1:3,13:15) = R_b_e * dtKf;
Phi_matrix(4:6,1:3) = -dtKf * crossProdMatrix(R_b_e * accMeasi);
Phi_matrix(4:6,4:6) = Phi_matrix(4:6,4:6) - 2 * Omega_ie * dtKf;
geocentric_radius = R_0 / sqrt(1 - (e * sin(latRad))^2) *...
    sqrt(cos(latRad)^2 + (1 - e^2)^2 * sin(latRad)^2); % from (2.137)
Phi_matrix(4:6,7:9) = -dtKf * 2 * Gravity_ECEF(pos) /...
    geocentric_radius * pos' / sqrt (pos' *...
    pos);
Phi_matrix(4:6,10:12) = R_b_e * dtKf;
Phi_matrix(7:9,4:6) = eye(3) * dtKf;
Phi_matrix(16,17) = dtKf;  
            
% 2. Determine approximate system noise covariance matrix using (14.82)
Q              = zeros(nState);
if PARAMS.stationaryMode
%     Q(1:3,1:3) = zeros(3);
%     Q(4:
    % don't add any process noise to the position, velocity, and attitude
elseif  (epoch-obj.epochLastInertialUpdate) >= dtKf 
    % We did not have IMU measurements available here
    Q(1:3,1:3)  = eye(3)*5^2;
    Q(4:6,4:6)  = eye(3)*5^2;
    Q(7:9,7:9)  = eye(3)*5^2;
else
    Q(1:3,1:3)     = eye(3) * PARAMS.gyro_noise_PSD * dtKf;
    Q(4:6,4:6)     = eye(3) * PARAMS.accel_noise_PSD * dtKf;
end
Q(10:12,10:12) = eye(3) * PARAMS.Q_ACC_BIAS * dtKf;
Q(13:15,13:15) = eye(3) * PARAMS.Q_W_BIAS * dtKf;
Q(16,16)       = PARAMS.Q_RXB.^2 * dtKf;
Q(17,17)       = PARAMS.Q_RXB.^2 * dtKf;

% 3. Propagate state estimates using (3.14) noting that only the clock
% states are non-zero due to closed-loop correction.
x_est_propagated = zeros(nState,1);
x_est_propagated(1:15,1) = 0;
x_est_propagated(16,1)   = 0;
x_est_propagated(17,1)   = 0;

% 4. Propagate state estimation error covariance matrix using (3.46)
cov_propagated = Phi_matrix * (cov + 0.5 * Q) *...
    Phi_matrix' + 0.5 * Q;
            
% Measurements are just position
pred_meas = zeros(6,1);
pred_meas(1:3) = pos-R_b_e*PARAMS.IMU_ARM;
pred_meas(4:6) = vel+R_b_e*crossProdMatrix(gyroMeasi)*PARAMS.IMU_ARM - ...
    Omega_ie*R_b_e*PARAMS.IMU_ARM;

% 5. Set-up measurement sensitivity matrix
H = zeros(6,nState);
H(1:3,7:9) = -eye(3);
H(1:3,4:6) = -eye(3);

% 6. Set-up measurement noise matrix
R = zeros(6);
R(1:3,1:3) = 1^2*eye(3);
R(4:6,4:6) = 0.2^2*eye(3);

% 7. Calculate Kalman gain using (3.21)
K = cov_propagated * H' * inv(H *cov_propagated * H' + R);

% 8. Measurement innovations
fullMeas = [posMeasi; velMeasi];
delta_z = fullMeas-pred_meas;

% 9. Update state estimates
x_est_new = x_est_propagated + K * delta_z;

% Update covariance
cov = (eye(nState) - K * H) * cov_propagated;

% CLOSED-LOOP CORRECTION

% Correct attitude, velocity, and position using (14.7-9)
R_b_e = (eye(3) - crossProdMatrix(x_est_new(1:3))) * R_b_e;
vel = vel - x_est_new(4:6);
pos = pos - x_est_new(7:9);

% Update IMU bias and GNSS receiver clock estimates
imuBiasStates = imuBiasStates + x_est_new(10:15);
% est_clock_new = x_est_new(16:17)';


% put updated values into object
obj.R_b_e = R_b_e;
obj.vel   = vel;
obj.pos   = pos;
obj.imuBiasStates = imuBiasStates;
obj.cov  = cov;

end






















