function outState = saveOutStatePpp(obj,outData,epoch,obs)
outState.epoch = epoch;

outState.pos = obj.pos;
outState.posApc = obj.posVelApc;

outState.covPos = obj.cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS);
outState.resids = obj.resids;
outState.residsInfo = [];
outState.measRemoved = obj.measRemoved;
outState.clockBias = obj.clockBias';

% Save tropospheric state information
outState.tropo = obj.state(obj.INDS_STATE.TROP);
outState.covTropo = obj.cov(obj.INDS_STATE.TROP,obj.INDS_STATE.TROP);

%
if ~isempty(obj.INDS_STATE.FLEX_STATES)
    % pull flex states
    
    flexStates = obj.state(obj.INDS_STATE.FLEX_STATES);
    flexStatesInfo = obj.INDS_STATE.FLEX_STATES_INFO;
    
    covDiag = diag(obj.cov);
    
    flexStatesCov = sqrt(covDiag(obj.INDS_STATE.FLEX_STATES));
    flexStatesEpochs = epoch*ones(size(flexStates));
    
    outState.flexStates.states = flexStates;
    outState.flexStates.info   = flexStatesInfo;
    outState.flexStates.std    = flexStatesCov;
    outState.flexStates.epochs = flexStatesEpochs;
else
    outState.flexStates = [];
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

end