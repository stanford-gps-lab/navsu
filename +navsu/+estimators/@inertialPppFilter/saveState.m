function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)


outState = [];

outState.epoch = epoch;
outState.pos   = obj.pos;
outState.resids = obj.resids;
outState.residsInfo = [];
outState.measRemoved = obj.measRemoved;

% don't include if this is an inertial measurement lol
measTypes = cellfun(@(x) getfield(x,'type'),obs);

if ismember(measTypes,navsu.internal.MeasEnum.IMU)
    return;
end

if isempty(outData) || isempty([outData(:).residsInfo])
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