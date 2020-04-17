function obs = measMask(obs,maskTable)

% maskTable = PARAMS.measMask;

maskVar = maskTable.Variables';

% for now, this assumes a fixed format for the measurement masking
rangeMask0 = reshape([maskVar(:,1) maskVar(:,2)]',size(maskVar,1)*2,1);
rangeMaskMat = repmat(rangeMask0,1,size(obs.range.obs,2));

doppMask0 = maskVar(1:3,3);
doppMaskMat = repmat(doppMask0,1,size(obs.doppler.obs,2));

%% Do the mask
obs.range.obs = obs.range.obs.*rangeMaskMat;
obs.doppler.obs = obs.doppler.obs.*doppMaskMat;


end