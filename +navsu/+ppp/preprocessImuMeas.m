function imuMeas = preprocessImuMeas(imuMeasRaw,varargin)

p = inputParser;

p.addParameter('epochStart',-Inf);
p.addParameter('epochEnd',Inf);
p.addParameter('downsampleFactor',1);

% parse the results
parse(p, varargin{:});
res = p.Results;
epochStart       = res.epochStart;       % Minimum time of observations
epochEnd         = res.epochEnd;         % Maximum time of observations
downsampleFactor = res.downsampleFactor; % Factor by which to downsample

if isempty(imuMeasRaw)
    imuMeas = [];
    return;
end

indsMeas = find(imuMeasRaw.epochs >= epochStart & imuMeasRaw.epochs < epochEnd);
indsMeas = indsMeas(1:downsampleFactor:end);

imuMeas.epochs      = imuMeasRaw.epochs(indsMeas,:);
imuMeas.headerTow   = imuMeasRaw.headerTow(indsMeas,:);
imuMeas.imuError    = imuMeasRaw.imuError(indsMeas,:);
imuMeas.imuType     = imuMeasRaw.imuType(indsMeas,:);
imuMeas.week        = imuMeasRaw.week(indsMeas,:);
imuMeas.tow         = imuMeasRaw.tow(indsMeas,:);
imuMeas.imuStatus   = imuMeasRaw.imuStatus(indsMeas,:);
imuMeas.acc         = imuMeasRaw.acc(indsMeas,:);
imuMeas.gyro        = imuMeasRaw.gyro(indsMeas,:);
imuMeas.type        = navsu.internal.MeasEnum.IMU;



end