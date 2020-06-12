function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)

% don't include if this is an inertial measurement lol
for jdx = 1:length(obs)
    measTypes(jdx) = obs{jdx}.type;
end
if ismember(measTypes,navsu.internal.MeasEnum.IMU)
    return;
end

outState = saveOutStatePpp(obj,outData,epoch,obs);

outState.posApc  = obj.posVelApc;
outState.pos = obj.pos;
outState.R_b_e = obj.R_b_e;
outState.vel   = obj.vel;

% Save accelerometer and gyro bias states
outState.imuBiasStates = obj.imuBiasStates;

% Save odometry state
outState.wheelScale = obj.state(obj.INDS_STATE.WHEELS);

if ~isempty(obj.lastGyroMeas)
    outState.gyro = obj.lastGyroMeas;
else
    outState.gyro = nan(3,1);
end

outData = [outData; outState];

end