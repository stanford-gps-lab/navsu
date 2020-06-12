function [predMeasi,Hi,Ri,measIdi,measi] = handlePositionMeas(obj,posMeas)


nState = size(obj.state,1);
if ~isempty(posMeas)
    if strcmp(posMeas.REFPOS,'REF')
        % IMU or otherwise reference location
        predMeasi = obj.pos;
    elseif strcmp(posMeas.REFPOS,'APC')
        predMeasi = obj.posVelApc;
    end
    
    Hi = zeros(3,nState);
    Hi(:,obj.INDS_STATE.POS) = -diag(ones(3,1));
    
    Ri = posMeas.cov;
    measIdi = posMeas.ID;
    measi = posMeas.obs';
else
    % Return empty stuff
    predMeasi = [];
    Hi = zeros(0,nState);
    Ri = [];
    measIdi = [];
    measi = [];
end

end