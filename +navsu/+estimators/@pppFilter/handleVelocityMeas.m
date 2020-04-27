function [predMeasi,Hi,Ri,measMati] = handleVelocityMeas(obj,velMeas)


nState = size(obj.state,1);
if ~isempty(velMeas)
    predMeasi = obj.vel;
    
    Hi = zeros(3,nState);
    Hi(:,obj.INDS_STATE.VEL) = -diag(ones(3,1));
    
    Ri = velMeas.cov;
    
    measMati = zeros(3,6);
    measMati(:,6) = 4;
    measMati(:,5) = velMeas.obs;
    
else
   % Return empty stuff
   predMeasi = [];
   Hi = zeros(0,nState);
   Ri = [];
   measMati = zeros(0,6);
    
end

end