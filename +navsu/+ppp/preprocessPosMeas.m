function posMeas = preprocessPosMeas(posRaw,varargin)

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

indsMeas = find(posRaw.epochs >= epochStart & posRaw.epochs < epochEnd);
indsMeas = indsMeas(1:downsampleFactor:end);

posMeas.epochs = posRaw.epochs(indsMeas,:);
posMeas.obs    = posRaw.obs(indsMeas,:);
posMeas.cov    = posRaw.cov(indsMeas,:,:);
posMeas.type   = posRaw.type;
posMeas.ID     = posRaw.ID;
posMeas.REFPOS = posRaw.REFPOS;

end