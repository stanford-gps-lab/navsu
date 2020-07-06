function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)

% don't include if this is an inertial measurement lol
measTypes = cellfun(@(x) getfield(x,'type'),obs);

if ismember(measTypes,navsu.internal.MeasEnum.IMU)
    return;
end

outState = saveOutStatePpp(obj,outData,epoch,obs);

% outState.pos   = obj.posVelApc;
outState.R_b_e = obj.R_b_e;
outState.vel   = obj.vel;

outState.imuBiasStates = obj.imuBiasStates;

if ~isempty(obj.lastGyroMeas)
    outState.gyro = obj.lastGyroMeas;
else
    outState.gyro = nan(3,1);
end

outData = [outData; outState];

end