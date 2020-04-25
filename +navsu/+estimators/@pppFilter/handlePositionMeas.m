function [predMeasi,Hi,Ri,measMati] = handlePositionMeas(obj,posMeas)


nState = size(obj.state,1);
if ~isempty(posMeas)
    
    
    predMeasi = obj.pos;
    
    Hi = zeros(3,nState);
    Hi(:,obj.INDS_STATE.POS) = -diag(ones(3,1));
    
    Ri = posMeas.cov;
    
    measMati = zeros(3,6);
    measMati(:,6) = 4;
    measMati(:,5) = posMeas.obs;
    
else
   % Return empty stuff
   predMeasi = [];
   Hi = zeros(0,nState);
   Ri = [];
   measMati = zeros(0,6);
    
end

end