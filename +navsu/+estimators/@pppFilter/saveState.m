function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)


outState = [];

outState.epoch = epoch;
if strcmp(obj.PARAMS.outputPos,'APC')
    outState.pos   = obj.pos;
elseif strcmp(obj.PARAMS.outputPos,'REF')
    % PPP filter does not update the attitude, so do so now
    if norm(obj.vel) > 0.5
        obj.R_b_e = navsu.ppp.posVel2Rbe(obj.pos,obj.vel);
    end
    outState.pos   = obj.pos-obj.R_b_e*obj.PARAMS.IMU_ARM*0;
    
end

outState.covPos = obj.cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS);
outState.resids = obj.resids;
outState.residsInfo = [];
outState.measRemoved = obj.measRemoved;

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