function wheels = preprocessWheels(wheelsRaw,varargin)

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

indsMeas = find(wheelsRaw.epochs >= epochStart & wheelsRaw.epochs < epochEnd);
indsMeas = indsMeas(1:downsampleFactor:end);

% Just looking for change in tick count- take diff for each 
dTickFrontLeft = [NaN; diff(wheelsRaw.SpeedFrontLeft)];
dTickFrontRight = [NaN; diff(wheelsRaw.SpeedFrontRight)];
dTickRearLeft = [NaN; diff(wheelsRaw.SpeedRearLeft)];
dTickRearRight = [NaN; diff(wheelsRaw.SpeedRearRight)];

wheels.epochs           = wheelsRaw.epochs(indsMeas,:);
wheels.headerTow        = wheelsRaw.headerTow(indsMeas,:);
wheels.SpeedFrontLeft   = dTickFrontLeft(indsMeas,:);
wheels.SpeedFrontRight  = dTickFrontRight(indsMeas,:);
wheels.SpeedRearLeft    = dTickRearLeft(indsMeas,:);
wheels.SpeedRearRight   = dTickRearRight(indsMeas,:);
wheels.SteeringAngle    = wheelsRaw.SteeringAngle(indsMeas,:);
wheels.VehicleSpeed     = wheelsRaw.VehicleSpeed(indsMeas,:);
wheels.Stationary       = wheelsRaw.Stationary(indsMeas,:);
wheels.TransmissionSetting  = wheelsRaw.TransmissionSetting(indsMeas,:);
wheels.ParkingBrakeStatus   = wheelsRaw.ParkingBrakeStatus(indsMeas,:);
wheels.headerTowNoRound  = wheelsRaw.headerTowNoRound(indsMeas,:);

wheels.type = navsu.internal.MeasEnum.Wheels;

end