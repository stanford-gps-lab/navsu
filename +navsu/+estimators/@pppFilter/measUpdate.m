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

llhi = navsu.geo.xyz2llh(pos');

cov_propagated  = obj.cov;
nState = size(cov_propagated,1);

epoch0 = epoch;

%% Pull off measurements
measMat = [];

% Wipe off measurements based on the desried measurement masking
obs = navsu.ppp.measMask(obs,PARAMS.measMask);

% Just using all PR measurements
indsMeasPr = find(obs.range.obs ~= 0  & obs.range.ind == 1 );
% indsMeasPr = [];
prnObsMat      = repmat(obs.PRN,size(obs.range.obs,1),1);
constIndObsMat = repmat(obs.constInds,size(obs.range.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnObsMat(indsMeasPr) constIndObsMat(indsMeasPr) ...
    obs.range.sig(indsMeasPr) obs.range.freqs(indsMeasPr) obs.range.obs(indsMeasPr) ...
    1*ones(size(indsMeasPr))]];

% Add carrier phase measurements
indsMeasPh = find(obs.range.obs ~= 0  & obs.range.ind == 2 );
% indsMeasPh = [];
prnObsMat      = repmat(obs.PRN,size(obs.range.obs,1),1);
constIndObsMat = repmat(obs.constInds,size(obs.range.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnObsMat(indsMeasPh) constIndObsMat(indsMeasPh) ...
    obs.range.sig(indsMeasPh) obs.range.freqs(indsMeasPh) obs.range.obs(indsMeasPh) ...
    2*ones(size(indsMeasPh))]];

% Add doppler measurements
indsMeasDop = find(obs.doppler.obs ~= 0);

prnDopMat      = repmat(obs.PRN,size(obs.doppler.obs,1),1);
constIndDopMat = repmat(obs.constInds,size(obs.doppler.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnDopMat(indsMeasDop) constIndDopMat(indsMeasDop) ...
    obs.doppler.sig(indsMeasDop) obs.doppler.freqs(indsMeasDop) ...
    obs.doppler.obs(indsMeasDop) 3*ones(size(indsMeasDop))]];

nMeas = size(measMat,1);

nSats = length(unique(measMat(:,1:2),'rows'));

if nMeas > 0
    %% Propagate orbit and clock for all measurements
    prnConstInds = sortrows(unique(measMat(:,1:2),'rows'),2);
    
    % transmission time for all satellites
    % rough estimate of travel time- just use geometric range
    [~,ib] = ismember(prnConstInds(:,2),obj.INDS_STATE.CLOCK_BIAS_CONSTS);
    tRx = epoch-obj.clockBias(ib)./c;
    [svPos,svVel] = corrData.propagate(prnConstInds(:,1),prnConstInds(:,2),tRx);
    
    % Might be missing some precise data- check and remove if so
    if any(isnan(svPos(:,1)))
        indsNan = find(isnan(svPos(:,1)));
        
        prnConstIndsNan = prnConstInds(indsNan,:);
        
        % Remove the associated measurements
        indsMeasRemove = find(ismember(measMat(:,1:2),prnConstIndsNan,'rows'));
        measMat(indsMeasRemove,:) = [];
        
        prnConstInds(indsNan,:) = [];
        [~,ib] = ismember(prnConstInds(:,2),obj.INDS_STATE.CLOCK_BIAS_CONSTS);
        tRx = epoch-obj.clockBias(ib)./c;
        
        svPos(indsNan,:) = [];
        
        nMeas = size(measMat,1);
    end
    
    satBias = -c*corrData.clock(prnConstInds(:,1),prnConstInds(:,2),tRx);
    gRangeSv = sqrt(sum((obj.pos'-svPos).^2,2));
    
    % Need rough estimate of the receiver clock bias in case of reset
    measMatPr = measMat(measMat(:,6) == 1,:);
    % Pull one pseudorange for each satellite
    [~,losIndPr] = ismember(measMatPr(:,1:2),prnConstInds,'rows');
    
    % bRxi = nanmedian(measi(measInfoi(:,4) == 1,1)- sqrt(sum((svPosRot-repmat(usrPos',size(svPosRot,1),1)).^2,2))-bSati(sIndsMap(measInfoi(:,4) == 1)));
    if isempty(measMatPr) || 1
        %         bRxi = obj.clockBias;
        bRxi = x_est_propagated(obj.INDS_STATE.CLOCK_BIAS,1);
    else
        bRxi = nanmedian(measMatPr(:,5)-gRangeSv(losIndPr)-satBias(losIndPr));
        obj.clockBias(:) = bRxi;
    end
    x_est_propagated(obj.INDS_STATE.CLOCK_BIAS,1)   = bRxi;
    
    epoch = epoch-obj.clockBias(ib)./c;
    
    satBias = -c*corrData.clock(prnConstInds(:,1),prnConstInds(:,2),epoch);
    [svPos,svVel] = corrData.propagate(prnConstInds(:,1),prnConstInds(:,2),epoch);
    gRangeSv = sqrt(sum((obj.pos'-svPos).^2,2));
    
    tTx = epoch-gRangeSv./c;
    
    [svPos,svVel] = corrData.propagate(prnConstInds(:,1),prnConstInds(:,2),tTx);
    travelTime = epoch-tTx;
    svPosRot = navsu.ppp.models.earthRotTravelCorr(travelTime,svPos);
    
    rxDrift   = obj.clockDrift(ib);
    rxBias    = obj.clockBias(ib);
    
    % elevation and azimuth for each LOS
    [el,az] = navsu.geo.pos2elaz(obj.pos',svPosRot);
    
    %% Various range effects
    % tropo delay for each LOS
    % trop = saastamoinen_model_SU(llhi(1), llhi(2), llhi(3), el*180/pi);
    % [trop0,m,tropDataSave] = tropo_error_correction_unb3(el,h,lat,doy);
    doy = navsu.time.jd2doy(navsu.time.epochs2jd(epoch));
    [trop,m,~] = navsu.ppp.models.tropDelay(el*180/pi,az*180/pi, llhi(:,3), llhi(:,1), llhi(:,2), doy, PARAMS, [],[],epoch);
    
    % TEC for each LOS
    if any(measMat(:,3) < 100 & measMat(:,6) < 3) %&& strcmp(PARAMS.states.ionoMode,'TEC')
        [~,~,tecSlant] = corrData.ionoDelay(epoch,llhi,'az',az,'el',el);
    else
        tecSlant = zeros(size(prnConstInds,1),1);
    end
    
    % solid tide adjustment to position
    [~,stRangeOffset] = navsu.ppp.models.solidTide(epoch(1),pos,'svPos',svPosRot);
    
    % relativistic corrections
    relClockCorr = 2/c^2.*sum(svPos.*svVel,2)*c;
    relRangeCorr = navsu.ppp.models.relRangeCorr(svPos',pos',PARAMS);
    
    % Carrier phase windup
    [~,ib] = ismember(prnConstInds,obj.phWind.PrnConstInd,'rows');
    phWind = navsu.ppp.models.carrierPhaseWindupGGM(epoch(1), repmat(pos',size(svPosRot,1)), svPosRot, obj.phWind.phaseOffset(ib));
    obj.phWind.phaseOffset(ib) = phWind; % need to update the phase windup object
    
    % Geometry matrix
    A = (svPosRot-pos')./sqrt(sum((pos'-svPosRot).^2,2));
    
    % geometric range
    gRange = sqrt(sum((svPosRot-pos').^2,2));
    dVel = vel'-svVel;
    
    %%
    [~,losInds] = ismember(measMat(:,[1 2]),prnConstInds,'rows');
    [~,indAmbStates] = ismember([measMat(:,[1 2 3]) ones(size(measMat,1),1)],obj.INDS_STATE.FLEX_STATES_INFO(:,[1 2 4 3]),'rows');
    [~,indIonos]   =ismember([measMat(:,1:2) 2*ones(size(measMat,1),1)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:3),'rows');
    [~,indGloDcbs] =ismember([measMat(:,1:2) 3*ones(size(measMat,1),1) measMat(:,3)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    [~,indMpCodes] =ismember([measMat(:,1:2) 4*ones(size(measMat,1),1) measMat(:,3)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    [~,indMpCarrs] =ismember([measMat(:,1:2) 5*ones(size(measMat,1),1) measMat(:,3)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    
    % build each measurement
    pred_meas = zeros(nMeas,1);
    H         = zeros(nMeas,nState);
    R         = zeros(nMeas,1);
    for idx = 1:nMeas
        measTypei = measMat(idx,6);% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
        prni      = measMat(idx,1);
        constIndi = measMat(idx,2);
        sigi      = measMat(idx,3);
        freqi     = measMat(idx,4);
        measi     = measMat(idx,5);
        
        losInd = losInds(idx);
        weighti = 1 ./ (sin(el(losInd)).^2);
        
        switch measTypei
            case 1
                % Code phase measurement
                [predMeasi,Hii,sigMeasi] = codeModel(obj,nState,sigi,freqi,tecSlant(losInd),x_est_propagated,...
                    constIndi,indGloDcbs(idx),indMpCodes(idx),m(losInd),gRange(losInd),satBias(losInd),rxBias(losInd),trop(losInd),stRangeOffset(losInd),...
                    relClockCorr(losInd),relRangeCorr(losInd),A(losInd,:));                
                
            case 2
                % Carrier phase measurement
                [predMeasi,Hii,sigMeasi] = carrierModel(obj,nState,sigi,freqi,...
                    tecSlant(losInd),x_est_propagated,m(losInd),indIonos(idx), ...
                    indMpCarrs(idx),indAmbStates(idx),phWind(losInd),...
                    gRange(losInd),satBias(losInd),rxBias(losInd),trop(losInd),...
                    stRangeOffset(losInd),relClockCorr(losInd),relRangeCorr(losInd),...
                    A(losInd,:),constIndi);
                
            case 3
                % Doppler measurement
                [predMeasi,Hii,sigMeasi] = obj.doppModel(nState,dVel(losInd,:),...
                    A(losInd,:),rxDrift(losInd),constIndi);
                
        end
        
        % Put computed predicted measurement, sensitivity, and
        % measurement sigma into their respective places.
        pred_meas(idx) = predMeasi;
        H(idx,:) = Hii;
        R(idx,idx) = weighti*sigMeasi^2;
    end
    
    % Remove measurements from satellites that are too low
    indsElLow = find(el(losInds)*180/pi < PARAMS.elMask);
    if ~isempty(indsElLow)
        measMatRemovedLow = measMat(indsElLow,:);
        
        H(indsElLow,:) = [];
        R(indsElLow,:) = [];
        R(:,indsElLow) = [];
        measMat(indsElLow,:) = [];
        pred_meas(indsElLow) = [];
    else
        measMatRemovedLow = zeros(0,size(measMat,2));
    end
else
    % No GNSS measurements
    % build each measurement
    pred_meas = zeros(0,1);
    H         = zeros(0,nState);
    R         = zeros(0,1);
    
end

if nMeas  == 0
    el = [];
    az = [];
    prnConstInds = zeros(0,2);
    measMatRemovedLow = zeros(0,2);
end

%% Pseudomeasurements
if PARAMS.measUse.noVertVel
    % No vertical velocity constraint (vehicle only moves forward)
    pseudoMeasi = 0;
    vertVeli = [0 0 1]*R_b_e'*velRx;
    
    ri = 1^2;
    Hi = zeros(1,size(H,2));
    Hi(1,obj.INDS_STATE.VEL) = -[0 0 1]*R_b_e';
    
    measMati = [0 0 0 0 pseudoMeasi 4];
    
    R(size(R,1)+1,size(R,1)+1) = ri;
    H = [H; Hi];
    measMat = [measMat; measMati];
    pred_meas = [pred_meas; vertVeli];
    
    pseudoMeasi = 0;
    vertVeli = [0 1 0]*R_b_e'*velRx;
    
    ri = 1^2;
    Hi = zeros(1,size(H,2));
    Hi(1,obj.INDS_STATE.VEL) = -[0 1 0]*R_b_e';
    
    measMati = [0 0 0 0 pseudoMeasi 4];
    
    R(size(R,1)+1,size(R,1)+1) = ri;
    H = [H; Hi];
    measMat = [measMat; measMati];
    pred_meas = [pred_meas; vertVeli];
end

nMeas = size(H,1);
%%
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


obj.resids.epoch = epoch0;
obj.resids.range = rangeResids;
obj.resids.doppler = doppResids;
obj.resids.el      = elFull;
obj.resids.az      = azFull;

obj.measRemoved.measRemove = measRemoveSave;
obj.measRemoved.epoch      = epochRemoveSave;

end