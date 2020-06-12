function [predMeasi,Hi,Ri,measIdi,measi] = handleVelocityMeas(obj,velMeas)


nState = size(obj.state,1);
if ~isempty(velMeas)
    if strcmp(velMeas.REFPOS,'REF')
        % IMU or otherwise reference location
        predMeasi = obj.vel;
    elseif strcmp(velMeas.REFPOS,'APC')
        [~,predMeasi] = obj.posVelApc;
    end
    
    Hi = zeros(3,nState);
    Hi(:,obj.INDS_STATE.VEL) = -diag(ones(3,1));
    
    Ri = velMeas.cov;
    measIdi = velMeas.ID;
    measi = velMeas.obs';
else
    % Return empty stuff
    predMeasi = [];
    Hi = zeros(0,nState);
    Ri = [];
    measIdi = [];
    measi = [];
end
end