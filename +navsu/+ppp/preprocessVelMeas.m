function velMeas = preprocessVelMeas(velRaw,varargin)

% this is mostly a wrapper for svPosFromProd
p = inputParser;

p.addParameter('epochStart',-Inf);
p.addParameter('epochEnd',Inf);
p.addParameter('downsampleFactor',1);

% parse the results
parse(p, varargin{:});
res = p.Results;
epochStart       = res.epochStart;       % Minimum time of observations
epochEnd         = res.epochEnd;         % Maximum time of observations
downsampleFactor = res.downsampleFactor; % FActor by which to downsample

if isempty(velRaw)
    velMeas = [];
    return;
end

indsMeas = find(velRaw.epochs >= epochStart & velRaw.epochs < epochEnd);
indsMeas = indsMeas(1:downsampleFactor:end);

velMeas.epochs = velRaw.epochs(indsMeas,:);
velMeas.obs    = velRaw.obs(indsMeas,:);
velMeas.cov    = velRaw.cov(indsMeas,:,:);
velMeas.type   = velRaw.type;
velMeas.ID     = velRaw.ID;
velMeas.REFPOS = velRaw.REFPOS;

end