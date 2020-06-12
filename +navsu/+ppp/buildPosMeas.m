function posMeas = buildPosMeas(epochs,pos,sig,idNum,varargin)

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

posMeas = [];

% Position input can be a single position.  If it is, need to make it match
% the size of 'epochs'
if size(pos,2) == 1
   pos = repmat(pos,1,length(epochs)); 
end

if all(size(sig) == [1 1])
   cov =  repmat(permute(diag(sig^2*ones(3,1)),[3 1 2]),length(epochs),1,1);
end

% Build position measurements
posMeas.epochs = epochs;
posMeas.obs    = pos;
posMeas.cov    = cov;
posMeas.type   = navsu.internal.MeasEnum.Position;
posMeas.ID     = navsu.internal.MeasIdPos(idNum*ones(3,1),...
    [navsu.internal.MeasEnum.PosX navsu.internal.MeasEnum.PosY navsu.internal.MeasEnum.PosZ]');
posMeas.REFPOS = REFPOS;

end