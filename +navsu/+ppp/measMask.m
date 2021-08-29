function obs = measMask(obs,maskTable)

% maskTable = PARAMS.measMask;

maskVar = maskTable.Variables';

% for now, this assumes a fixed format for the measurement masking
rangeMask = reshape(maskVar(:, 1:2)', size(maskVar,1)*2, 1);

doppMask = maskVar(1:3, 3); % single freq only

%% Do the mask
obs.range.obs = obs.range.obs .* rangeMask;
obs.doppler.obs = obs.doppler.obs .* doppMask;


end