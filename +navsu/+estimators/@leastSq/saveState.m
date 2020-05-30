function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)


outState = [];

outState.epoch = epoch;
outState.pos   = obj.pos; 
outState.covPos = obj.cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS);

outData = [outData; outState];

end