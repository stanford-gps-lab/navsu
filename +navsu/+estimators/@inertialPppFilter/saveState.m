function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)


outState = [];


% don't include if this is an inertial measurement lol
measTypes = cellfun(@(x) getfield(x,'type'),obs);

if ismember(measTypes,navsu.internal.MeasEnum.IMU)
    return;
end

outState.epoch = epoch;
outState.pos   = obj.posVelApc;
outState.resids = obj.resids;
outState.residsInfo = [];
outState.measRemoved = obj.measRemoved;
outState.covPos = obj.cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS);
outState.R_b_e = obj.R_b_e;
outState.vel   = obj.vel;

if ~isempty(obj.lastGyroMeas)
    outState.gyro = obj.lastGyroMeas;
else
    outState.gyro = nan(3,1);
end

if isempty(outData) || isempty([outData(:).residsInfo]) && ~ ismember(measTypes,navsu.internal.MeasEnum.IMU)
    gnssMeas = [];
    for idx = 1:length(obs)
        obsi = obs{idx};
        switch obsi.type
            case navsu.internal.MeasEnum.GNSS
                gnssMeas = obsi;
        end
    end
    
    if ~isempty(gnssMeas)
        
        % this is the first one- include some more information
        outState.residsInfo.rangeInfo = gnssMeas.range;
        outState.residsInfo.rangeInfo.obs = [];
        outState.residsInfo.rangeInfo.lockTime = [];
        
        outState.residsInfo.dopplerInfo = gnssMeas.doppler;
        outState.residsInfo.dopplerInfo.obs = [];
    end
end

outData = [outData; outState];

end