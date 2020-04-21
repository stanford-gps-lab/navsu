function measUpdate(obj,epoch,obs,corrData,measRemovedSlip)

%% Pull a few things out of the filter object

% Run parameters
PARAMS = obj.PARAMS;

% Speed of light
c = navsu.constants.c;

% State- has been propagated in the time update
x_est_propagated = obj.state;

% Receiver position
pos = obj.pos;

% Receiver velocity
vel = obj.vel;

% 
cov_propagated  = obj.cov;
nState = size(cov_propagated,1);

epoch0 = epoch;

%% GNSS Measurements
[pred_meas,H,R,el,az,prnConstInds,measMatRemovedLow,measMat] = handleGnssMeas(obj,epoch0,obs,corrData);

%% Pseudomeasurements
if PARAMS.measUse.noVertVel
   [predMeasi,Hi,Ri,measMati] = handleVehicleConstraintPseudomeas(obj);

    % Add everything
    pred_meas = [pred_meas; predMeasi];
    H         = [H; Hi];
    measMat   = [measMat; measMati];
    
    R2 = zeros(size(R,1)+2,size(R,1)+2);
    R2(1:size(R,1),1:size(R,1)) = R;
    R2((end-1):end,(end-1):end) = Ri;
    R = R2;
end

nMeas = size(H,1);
%% Do the measurement update
if nMeas > 0
    % Measurement were available- do the update.
    [H,delta_z,residsPost,K,measMat,~,measMatRemoved,R] = ...
        navsu.ppp.measUpdateExclude(H,cov_propagated,R,measMat,pred_meas,PARAMS);
    
    % 9. Update state estimates
    x_est_new = x_est_propagated + K * delta_z;
    
    % Update covariance
    cov = (eye(nState) - K * H) * cov_propagated;
    
else
    % No measurement update
    x_est_new = x_est_propagated;
    
    cov = cov_propagated;
    
    measMatRemoved = zeros(0,6);
    residsPost = [];
end

%% Update the position and velocity values
vel = vel - x_est_new(obj.INDS_STATE.VEL);
pos = pos - x_est_new(obj.INDS_STATE.POS);

obj.clockBias  = x_est_new(obj.INDS_STATE.CLOCK_BIAS);
obj.clockDrift = x_est_new(obj.INDS_STATE.CLOCK_DRIFT);

% put updated values into object
obj.vel   = vel;
obj.pos   = pos;
obj.cov  = cov;
obj.state = x_est_new;
obj.posPrevTc = pos;

% Deal with resets if any of the removed measurements were carrier phases
if ~isempty(measMatRemoved)
    for jdx = 1:size(measMatRemoved,1)
        if measMatRemoved(jdx,6) == 2 %|| measMatRemoved(jdx,6) == 1
            % carrier phase- reset ambiguity by just removing the state
            obj.removeFlexState([measMatRemoved(jdx,[1 2]) 1 measMatRemoved(jdx,3)] );
        end
    end
end

%% Save things for output
obj.allSatsSeen = sortrows(unique([measMat(:,1:2); obj.allSatsSeen],'rows'),2);

[~,indsSave] = ismember(measMat(measMat(:,end) == 1 | ...
    measMat(:,end) == 2,[1 2 3 6]),[obs.range.PRN(:) obs.range.constInds(:) ...
    obs.range.sig(:) obs.range.ind(:)],'rows');
rangeResids = nan(size(obs.range.obs));
rangeResids(indsSave) = residsPost(measMat(:,6) == 1 | measMat(:,6) == 2);

% Save the doppler residuals
[~,indsSave] = ismember(measMat(measMat(:,end) == 3 ,[1 2 3]),...
    [obs.doppler.PRN(:) obs.doppler.constInds(:) ...
    obs.doppler.sig(:) ],'rows');
doppResids = nan(size(obs.doppler.obs));
doppResids(indsSave) = residsPost(measMat(:,6) == 3);

[~,indsEl] = ismember(prnConstInds,[obs.PRN' obs.constInds'],'rows');

elFull = nan(size(el,1),1);
elFull(indsEl) = el;
azFull = nan(size(el,1),1);
azFull(indsEl) = az;

% Actually only keep the prn, const, and reason for elevation removals
measLow = unique(measMatRemovedLow(:,[1 2]),'rows');

% PRN | CONST | SIG  | MEAS TYPE (1=CODE,2=CARR,3=DOP) | REMOVAL REASON (1=ELEVATION,2=RESIDUALS, 3 = SLIP)
measRemoveSave = [measLow(:,1:2) 0*ones(size(measLow,1),2) 1*ones(size(measLow,1),1);
    measMatRemoved(:,[1:3 6]) 2*ones(size(measMatRemoved,1),1);
    measRemovedSlip(:,[1 2 4]) 2*ones(size(measRemovedSlip,1),1) 3*ones(size(measRemovedSlip,1),1)];

epochRemoveSave = epoch0*ones(size(measRemoveSave,1),1);

%% Put everything back in the filter object.
obj.resids.epoch = epoch0;
obj.resids.range = rangeResids;
obj.resids.doppler = doppResids;
obj.resids.el      = elFull;
obj.resids.az      = azFull;
obj.measRemoved.measRemove = measRemoveSave;
obj.measRemoved.epoch      = epochRemoveSave;

end