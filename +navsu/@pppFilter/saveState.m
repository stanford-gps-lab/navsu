function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)


outState = [];

outState.epoch = epoch;
outState.pos   = obj.pos; 
outState.resids = obj.resids;
outState.residsInfo = [];

if isempty(outData)
    % this is the first one- include some more information
    outState.residsInfo.rangeInfo = obs.range;
    outState.residsInfo.rangeInfo.obs = [];
    outState.residsInfo.rangeInfo.lockTime = [];
    
    outState.residsInfo.dopplerInfo = obs.doppler;
    outState.residsInfo.dopplerInfo.obs = [];
    
end

outData = [outData; outState];


end