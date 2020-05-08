function [predMeasi,Hi,Ri,measIdi,measi] = handleVelocityMeas(obj,velMeas)


nState = size(obj.state,1);
if ~isempty(velMeas)
    predMeasi = obj.vel;
    
    Hi = zeros(3,nState);
    Hi(:,obj.INDS_STATE.VEL) = -diag(ones(3,1));
    
    Ri = velMeas.cov;
    
    measIdi = velMeas.ID;
    measi = velMeas.obs;
    
else
    % Return empty stuff
    predMeasi = [];
    Hi = zeros(0,nState);
    Ri = [];
    measIdi = [];
    measi = [];

end

end