function obs = checkMeas(obj,obs0)

% Go through each obs and only keep what's usable for this filter

% TO DO: also exclude measurements that are turned off in the filter
% parameters


% Check for IMU measurements- these are actually REQUIRED for the IMU
% filter.
measTypes = cellfun(@(x) getfield(x,'type'),obs0);
if ~ismember(navsu.internal.MeasEnum.IMU,measTypes)
    error('IMU measurements required for this filter type')
end

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
           obs = [obs obs0(idx)];
    end    
end
warning('on','backtrace')
