function velMeas = buildVelMeas(epochs,vel,sig,idNum)

velMeas = [];

% Position input can be a single position.  If it is, need to make it match
% the size of 'epochs'
if size(vel,2) == 1
   vel = repmat(vel,1,length(epochs)); 
end

if all(size(sig) == [1 1])
   cov =  repmat(diag(sig^2*ones(3,1)),1,1,length(epochs));
end

% Build position measurements
velMeas.epochs = epochs;
velMeas.obs    = vel;
velMeas.cov    = cov;
velMeas.type   = navsu.internal.MeasEnum.Velocity;
velMeas.ID     = navsu.internal.MeasIdVel(idNum*ones(3,1),...
    [navsu.internal.MeasEnum.VelX navsu.internal.MeasEnum.VelY navsu.internal.MeasEnum.VelZ]');

end