function velMeas = buildVelMeas(epochs,vel,sig,idNum,varargin)

% epochs and position must have time indexed in the rows. (i.e. column
% vectors) 

p = inputParser;

expectedRefPos = {'REF' 'APC'};
p.addParameter('REFPOS','REF', @(x) any(validatestring(x,expectedRefPos))); 

% parse the results
parse(p, varargin{:});
res = p.Results;
% position measurement can refer to the 'reference position' i.e. the IMU
% location or the antenna phase center 'REF' or 'APC' inputs respectively
REFPOS          = res.REFPOS;  

velMeas = [];

% Position input can be a single velocity (if the antenna is stationary). 
% If it is, need to make it match the size of 'epochs'
if size(vel,2) == 1
   vel = repmat(vel,1,length(epochs)); 
end

if all(size(sig) == [1 1])
   cov =  repmat(permute(diag(sig^2*ones(3,1)),[3 1 2]),length(epochs),1,1);
end

% Build position measurements
velMeas.epochs = epochs;
velMeas.obs    = vel;
velMeas.cov    = cov;
velMeas.type   = navsu.internal.MeasEnum.Velocity;
velMeas.ID     = navsu.internal.MeasIdVel(idNum*ones(3,1),...
    [navsu.internal.MeasEnum.VelX navsu.internal.MeasEnum.VelY navsu.internal.MeasEnum.VelZ]');
velMeas.REFPOS = REFPOS;


end