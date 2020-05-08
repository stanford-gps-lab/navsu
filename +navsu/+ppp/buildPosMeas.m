function posMeas = buildPosMeas(epochs,pos,sig,idNum)

posMeas = [];

% Position input can be a single position.  If it is, need to make it match
% the size of 'epochs'
if size(pos,2) == 1
   pos = repmat(pos,1,length(epochs)); 
end

if all(size(sig) == [1 1])
   cov =  repmat(diag(sig^2*ones(3,1)),1,1,length(epochs));
end

% Build position measurements
posMeas.epochs = epochs;
posMeas.obs    = pos;
posMeas.cov    = cov;
posMeas.type   = navsu.internal.MeasEnum.Position;
posMeas.ID     = navsu.internal.MeasIdPos(idNum*ones(3,1),...
    [navsu.internal.MeasEnum.PosX navsu.internal.MeasEnum.PosY navsu.internal.MeasEnum.PosZ]');

end