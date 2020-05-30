function obs = checkMeas(obj,obs0)

% Go through each obs and only keep what's usable for this filter

% For the standard PPP filter, exclude IMU measurements.

% TO DO: also exclude measurements that are turned off in the filter
% parameters

warning('off','backtrace')

obs = [];
for idx = 1:length(obs0)
    switch obs0{idx}.type
        case navsu.internal.MeasEnum.GNSS
            obs = [obs obs0(idx)];
            
        case navsu.internal.MeasEnum.Position
            obs = [obs obs0(idx)];
            
        case navsu.internal.MeasEnum.Velocity
            obs = [obs obs0(idx)];
            
        case navsu.internal.MeasEnum.IMU
            warning(['Excluding IMU measurements- not available in this' ...
                ' standard PPP filter- maybe try navsu.estimators.inertialPppFilter']);
            continue;
    end    
end
warning('on','backtrace')
