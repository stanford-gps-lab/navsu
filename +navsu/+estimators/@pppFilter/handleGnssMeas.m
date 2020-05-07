function  [pred_meas,H,R,el,az,prnConstInds,measMatRemovedLow,measMat] = handleGnssMeas(obj,epoch,obs,corrData)


measMat = [];

PARAMS = obj.PARAMS;

c = navsu.constants.c;

x_est_propagated = obj.state;

% Antenna phase center position and velocity
[pos,vel] = obj.posVelApc;

% pos = obj.pos;

% Receiver velocity
% vel = obj.vel;

llhi = navsu.geo.xyz2llh(pos');

nState = size(x_est_propagated,1);

%% Pull off measurements
% Wipe off measurements based on the desried measurement masking
obs = navsu.ppp.measMask(obs,PARAMS.measMask);

% Just using all PR measurements
indsMeasPr = find(obs.range.obs ~= 0  & obs.range.ind == 1 );
% indsMeasPr = [];
prnObsMat      = repmat(obs.PRN,size(obs.range.obs,1),1);
constIndObsMat = repmat(obs.constInds,size(obs.range.obs,1),1);

idList = [];
measList = [];

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnObsMat(indsMeasPr) constIndObsMat(indsMeasPr) ...
    obs.range.sig(indsMeasPr) obs.range.freqs(indsMeasPr) obs.range.obs(indsMeasPr) ...
    1*ones(size(indsMeasPr))]];
idList = [idList; obs.range.ID(indsMeasPr)];
measList = [measList; obs.range.obs(indsMeasPr)];

% Add carrier phase measurements
indsMeasPh = find(obs.range.obs ~= 0  & obs.range.ind == 2 );
% indsMeasPh = [];
prnObsMat      = repmat(obs.PRN,size(obs.range.obs,1),1);
constIndObsMat = repmat(obs.constInds,size(obs.range.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnObsMat(indsMeasPh) constIndObsMat(indsMeasPh) ...
    obs.range.sig(indsMeasPh) obs.range.freqs(indsMeasPh) obs.range.obs(indsMeasPh) ...
    2*ones(size(indsMeasPh))]];
idList = [idList; obs.range.ID(indsMeasPh)];
measList = [measList; obs.range.obs(indsMeasPh)];


% Add doppler measurements
indsMeasDop = find(obs.doppler.obs ~= 0);

prnDopMat      = repmat(obs.PRN,size(obs.doppler.obs,1),1);
constIndDopMat = repmat(obs.constInds,size(obs.doppler.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnDopMat(indsMeasDop) constIndDopMat(indsMeasDop) ...
    obs.doppler.sig(indsMeasDop) obs.doppler.freqs(indsMeasDop) ...
    obs.doppler.obs(indsMeasDop) 3*ones(size(indsMeasDop))]];

idList = [idList; obs.doppler.ID(indsMeasDop)];
measList = [measList; obs.doppler.obs(indsMeasDop)];


nMeas = size(measMat,1);

if nMeas > 0
    %% Propagate orbit and clock for all measurements
%     prnConstInds = sortrows(unique(measMat(:,1:2),'rows'),2);
    prnConstInds = sortrows(unique([[idList.prn]' [idList.const]'],'rows'),2); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        indsMeasRemove = find(ismember([[idList.prn]' [idList.const]'],prnConstIndsNan,'rows'));
        measMat(indsMeasRemove,:) = [];
        idList(indsMeasRemove) = [];
        measList(indsMeasRemove) = [];
        
        prnConstInds(indsNan,:) = [];
        [~,ib] = ismember(prnConstInds(:,2),obj.INDS_STATE.CLOCK_BIAS_CONSTS);
        tRx = epoch-obj.clockBias(ib)./c;
        
        svPos(indsNan,:) = [];
        
%         nMeas = size(measMat,1);
        nMeas = length(idList);
    end
    
    satBias = -c*corrData.clock(prnConstInds(:,1),prnConstInds(:,2),tRx);
    gRangeSv = sqrt(sum((obj.pos'-svPos).^2,2));
    
    % Need rough estimate of the receiver clock bias in case of reset
    measMatPr0 = measMat(measMat(:,6) == 1,:);
    indPr = find([idList.subtype]' == 1);
    measPr = measList(indPr,:);
    idPr   = idList(indPr);
    % Pull one pseudorange for each satellite
    [~,losIndPr] = ismember([[idPr.prn]' [idPr.const]'],prnConstInds,'rows');
    
    if isempty(measPr) || 1
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




end