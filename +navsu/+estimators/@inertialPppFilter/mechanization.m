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

end




