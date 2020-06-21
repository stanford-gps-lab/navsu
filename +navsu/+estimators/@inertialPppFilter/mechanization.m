function mechanization(obj,epoch,obs)

% Pull off just the IMU meas
imuMeas = navsu.ppp.pullMeasFromList(obs,navsu.internal.MeasEnum.IMU);

if ~isempty(imuMeas)
    accMeas  = imuMeas.acc' - obj.imuBiasStates(1:3);
    gyroMeas = imuMeas.gyro' - obj.imuBiasStates(4:6);
    
    obj.lastAccMeas = accMeas;
    obj.lastGyroMeas = gyroMeas;
    
else
    % Use our last acceleration and gyro measurements
    accMeas = obj.lastAccMeas;
    gyroMeas = obj.lastGyroMeas;
end


posOld = obj.pos;
velOld = obj.vel;
attOld = obj.R_b_e;

% epoch =

omegaIe = navsu.constants.omegaEarth; % earth rotation rate in rad/s

% time increment in seconds
dt = min([epoch-obj.epochLastInertialUpdate epoch-obj.epochLastGnssUpdate]);
obj.epochLastInertialUpdate = epoch;

% This mechanization is from Groves :)

% attitude update
alphaIe = omegaIe*dt;
% earth rotation matrix during the interval
cEarth = [cos(alphaIe), sin(alphaIe), 0;...
    -sin(alphaIe), cos(alphaIe), 0;...
    0,             0,  1];


% Compute attitude increment
alphaIbB = gyroMeas*dt;
magAlpha = norm(alphaIbB);
alphaIbBMat = navsu.geo.crossProdMatrix(alphaIbB);

% Obtain coordinate transformation matrix from the new attitude w.r.t. an
% inertial frame to the old using Rodrigues' formula, (5.73)
if magAlpha>1.E-8
    CNewOld = eye(3) + sin(magAlpha) / magAlpha * alphaIbBMat +...
        (1 - cos(magAlpha)) / magAlpha^2 * alphaIbBMat * alphaIbBMat;
else
    CNewOld = eye(3) + alphaIbBMat;
end


% Update attitude using (5.75)
att = cEarth * attOld * CNewOld;

% SPECIFIC FORCE FRAME TRANSFORMATION
% Calculate the average body-to-ECEF-frame coordinate transformation
% matrix over the update interval using (5.84) and (5.85)
if magAlpha>1.E-8
    aveCbe = attOld * (eye(3) + (1 - cos(magAlpha)) / magAlpha^2 ...
        * alphaIbBMat + (1 - sin(magAlpha) / magAlpha) / magAlpha^2 ...
        * alphaIbBMat * alphaIbBMat) - 0.5 * navsu.geo.crossProdMatrix([0;0;alphaIe])...
        * attOld;
else
    aveCbe = attOld - 0.5 * navsu.geo.crossProdMatrix([0;0;alphaIe]) *...
        attOld;
end %if mag_alpha


% Transform specific force to ECEF-frame resolving axes using (5.85)
fIbe = aveCbe * accMeas;

% UPDATE VELOCITY
% From (5.36),
vel = velOld + dt * (fIbe + navsu.geo.gravityEcef(posOld) -...
    2 * navsu.geo.crossProdMatrix([0;0;omegaIe]) * velOld);

% UPDATE CARTESIAN POSITION
% From (5.38),
pos = posOld + (vel + velOld) * 0.5 * dt;

% Update our object
obj.pos = pos;
obj.vel = vel;
obj.R_b_e = att;

velMean = (velOld+vel)/2;

% pitching teh car frame?
% roll pitch yaw
R2 = navsu.geo.euler2dcm123([0 0.8 0]*pi/180);

% rear left wheel lever arm
l_rl = obj.PARAMS.ARM_REF_REAR_AXLE + ...
    [0 -1/2*obj.PARAMS.WHEEL_TRACK -1/2*obj.PARAMS.WHEEL_DIAM]';

% rear right wheel lever arm
l_rr = obj.PARAMS.ARM_REF_REAR_AXLE + ...
    [0 1/2*obj.PARAMS.WHEEL_TRACK -1/2*obj.PARAMS.WHEEL_DIAM]';

% If using wheel odometry, integrate some things here:
if obj.PARAMS.states.wheels
    % rear left wheel velocity estimate increment    
    vel_rl_dt = [-1 0 0]*R2*(att'*velMean-navsu.geo.crossProdMatrix(gyroMeas)*l_rl)*dt;

    % rear right wheel velocity estimate increment
    vel_rr_dt = [-1 0 0]*R2*(att'*velMean-navsu.geo.crossProdMatrix(gyroMeas)*l_rr)*dt;
    
    % attitude sensitivity increment
    H11_dt = [-1 0 0]*R2*att'*navsu.geo.crossProdMatrix(velMean)*dt;
    
    % velocity sensitivity increment
    H12_dt = [-1 0 0]*R2*att'*dt;
        
    % Add all of these
    obj.wheelInfo.vrl_int = obj.wheelInfo.vrl_int + vel_rl_dt;
    obj.wheelInfo.vrr_int = obj.wheelInfo.vrr_int + vel_rr_dt;
    obj.wheelInfo.H11_int = obj.wheelInfo.H11_int + H11_dt;
    obj.wheelInfo.H12_int = obj.wheelInfo.H12_int + H12_dt;
    obj.wheelInfo.dt_int = obj.wheelInfo.dt_int+dt;

end

if obj.PARAMS.measUse.noVertVel
     % vertical velocity estimate at rear right wheel
    vel_up_dt = [0 0 1]*R2*(att'*velMean-navsu.geo.crossProdMatrix(gyroMeas)*l_rr)*dt;
    
    % vertical velocity sensitivities
    Hv1_dt = [0 0 1]*R2*att'*navsu.geo.crossProdMatrix(velMean)*dt;
    Hv2_dt = [0 0 1]*R2*att'*dt;
    
    % cross track velocity estimate at rear right wheel
    vel_cr_dt = [0 1 0]*R2*(att'*velMean-navsu.geo.crossProdMatrix(gyroMeas)*l_rr)*dt;
    
    % cross track velocity sensitivities
    Hc1_dt = [0 1 0]*R2*att'*navsu.geo.crossProdMatrix(velMean)*dt;
    Hc2_dt = [0 1 0]*R2*att'*dt;
    
    obj.wheelInfo.vup_int = obj.wheelInfo.vup_int + vel_up_dt;
    obj.wheelInfo.Hv1_int = obj.wheelInfo.Hv1_int + Hv1_dt;
    obj.wheelInfo.Hv2_int = obj.wheelInfo.Hv2_int + Hv2_dt;
    
    obj.wheelInfo.vcr_int = obj.wheelInfo.vcr_int + vel_cr_dt;
    obj.wheelInfo.Hc1_int = obj.wheelInfo.Hc1_int + Hc1_dt;
    obj.wheelInfo.Hc2_int = obj.wheelInfo.Hc2_int + Hc2_dt;
    
end


end




