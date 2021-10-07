function wheels = preprocessWheels(wheelsRaw, varargin)

% parse optional inputs
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

wheels.epochs           = wheelsRaw.epochs(indsMeas,:);
wheels.headerTow        = wheelsRaw.headerTow(indsMeas,:);
wheels.SpeedFrontLeft   = wheelsRaw.SpeedFrontLeft(indsMeas,:);
wheels.SpeedFrontRight  = wheelsRaw.SpeedFrontRight(indsMeas,:);
wheels.SpeedRearLeft    = wheelsRaw.SpeedRearLeft(indsMeas,:);
wheels.SpeedRearRight   = wheelsRaw.SpeedRearRight(indsMeas,:);
wheels.SteeringAngle    = wheelsRaw.SteeringAngle(indsMeas,:);
wheels.VehicleSpeed     = wheelsRaw.VehicleSpeed(indsMeas,:);
wheels.Stationary       = wheelsRaw.Stationary(indsMeas,:);
wheels.TransmissionSetting  = wheelsRaw.TransmissionSetting(indsMeas,:);
wheels.ParkingBrakeStatus   = wheelsRaw.ParkingBrakeStatus(indsMeas,:);

wheels.type = navsu.internal.MeasEnum.Wheels;

end