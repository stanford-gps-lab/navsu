function [predMeasi,Hi,Ri,measIdi,measi] = handlePositionMeas(obj,posMeas)


nState = size(obj.state,1);
if ~isempty(posMeas)
    predMeasi = obj.pos;
    
    Hi = zeros(3,nState);
    Hi(:,obj.INDS_STATE.POS) = -diag(ones(3,1));
    
    Ri = posMeas.cov;
    measIdi = posMeas.ID;
    measi = posMeas.obs;
else
    % Return empty stuff
    predMeasi = [];
    Hi = zeros(0,nState);
    Ri = [];
    measIdi = [];
    measi = [];
end

end