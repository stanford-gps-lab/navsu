function  [predMeas,H,R,el,az,prnConstInds,idList,measList,measIdRemovedLow] ...
    = handleGnssMeas(obj,epoch,obs,corrData,varargin)

p = inputParser;

p.addParameter('SimpleModel',false);

% parse the results
parse(p, varargin{:});
res        = p.Results;
SimpleModel= res.SimpleModel;

%%

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
% Wipe off measurements based on the desired measurement masking
obs = navsu.ppp.measMask(obs,PARAMS.measMask);

% Just using all PR measurements
indsMeasPr = find(obs.range.obs ~= 0  & obs.range.ind == 1 );

idList = [];
measList = [];

idList = [idList; obs.range.ID(indsMeasPr)];
measList = [measList; obs.range.obs(indsMeasPr)];

% Add carrier phase measurements
indsMeasPh = find(obs.range.obs ~= 0  & obs.range.ind == 2 );
idList = [idList; obs.range.ID(indsMeasPh)];
measList = [measList; obs.range.obs(indsMeasPh)];

% Add doppler measurements
indsMeasDop = find(obs.doppler.obs ~= 0);
idList = [idList; obs.doppler.ID(indsMeasDop)];
measList = [measList; obs.doppler.obs(indsMeasDop)];

nMeas = length(idList);

if nMeas > 0
    %% Propagate orbit and clock for all measurements
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
        idList(indsMeasRemove) = [];
        measList(indsMeasRemove) = [];
        
        prnConstInds(indsNan,:) = [];
        [~,ib] = ismember(prnConstInds(:,2),obj.INDS_STATE.CLOCK_BIAS_CONSTS);
        tRx = epoch-obj.clockBias(ib)./c;
        
        svPos(indsNan,:) = [];
        
        nMeas = length(idList);
    end
    
    satBias = -c*corrData.clock(prnConstInds(:,1),prnConstInds(:,2),tRx);
    gRangeSv = sqrt(sum((obj.pos'-svPos).^2,2));
    
    % Need rough estimate of the receiver clock bias in case of reset
    indPr = find([idList.subtype]' == navsu.internal.MeasEnum.Code);
    measPr = measList(indPr,:);
    idPr   = idList(indPr);
    % Pull one pseudorange for each satellite
    [~,losIndPr] = ismember([[idPr.prn]' [idPr.const]'],prnConstInds,'rows');
    
    if isempty(measPr) || 1
        bRxi = x_est_propagated(obj.INDS_STATE.CLOCK_BIAS,1);
    else
%         bRxi = nanmedian(measMatPr(:,5)-gRangeSv(losIndPr)-satBias(losIndPr));
%         obj.clockBias(:) = bRxi;
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
    if norm(pos) < 1e3
        trop = zeros(size(el));
        m = zeros(size(el));
        tecSlant = zeros(size(el));
    else
        doy = navsu.time.jd2doy(navsu.time.epochs2jd(epoch));
        [trop,m,~] = navsu.ppp.models.tropDelay(el*180/pi,az*180/pi, llhi(:,3), llhi(:,1), llhi(:,2), doy, PARAMS, [],[],epoch);
        
        % TEC for each LOS
        if any([idList.freq] < 100 & (([idList.subtype] == navsu.internal.MeasEnum.Code) | ([idList.subtype] == navsu.internal.MeasEnum.Carrier))) %&& strcmp(PARAMS.states.ionoMode,'TEC')
            [~,~,tecSlant] = corrData.ionoDelay(epoch,llhi,'az',az,'el',el);
        else
            tecSlant = zeros(size(prnConstInds,1),1);
        end
    end
    
    % solid tide adjustment to position
    if norm(pos) < 1e3
        stRangeOffset = zeros(size(svPosRot,1),1);
    else
        [~,stRangeOffset] = navsu.ppp.models.solidTide(epoch(1),pos,'svPos',svPosRot);
    end
    
    % relativistic corrections
    relClockCorr = 2/c^2.*sum(svPos.*svVel,2)*c;
    if norm(pos) < 1e3
        relRangeCorr = zeros(size(svPos,1),1);
    else 
        relRangeCorr = navsu.ppp.models.relRangeCorr(svPos',pos',PARAMS);
    end
    
    % Carrier phase windup
    [~,ib] = ismember(prnConstInds,obj.phWind.PrnConstInd,'rows');
    if all(ib)
        phWind = navsu.ppp.models.carrierPhaseWindupGGM(epoch(1), repmat(pos',size(svPosRot,1)), svPosRot, obj.phWind.phaseOffset(ib));
        obj.phWind.phaseOffset(ib) = phWind; % need to update the phase windup object
    end
    
    % Geometry matrix
    A = (svPosRot-pos')./sqrt(sum((pos'-svPosRot).^2,2));
    
    % geometric range
    gRange = sqrt(sum((svPosRot-pos').^2,2));
    dVel = vel'-svVel;
    
    %%
    
    [~,losInds] = ismember([[idList.prn]' [idList.const]'],prnConstInds,'rows');
    [~,indAmbStates] = ismember([[[idList.prn]' [idList.const]' [idList.freq]'] ones(length(idList),1)],obj.INDS_STATE.FLEX_STATES_INFO(:,[1 2 4 3]),'rows');
    [~,indIonos]   =ismember([[[idList.prn]' [idList.const]'] 2*ones(length(idList),1)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:3),'rows');
    [~,indGloDcbs] =ismember([[[idList.prn]' [idList.const]'] 3*ones(length(idList),1) [idList.freq]'],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    [~,indMpCodes] =ismember([[[idList.prn]' [idList.const]'] 4*ones(length(idList),1) [idList.freq]'],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    [~,indMpCarrs] =ismember([[[idList.prn]' [idList.const]'] 5*ones(length(idList),1) [idList.freq]'],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    
    
    % build each measurement
    predMeas = zeros(nMeas,1);
    H         = zeros(nMeas,nState);
    R         = zeros(nMeas,1);
    for idx = 1:nMeas
        measTypei = idList(idx).subtype;
        prni      = idList(idx).prn;
        constIndi = idList(idx).const;
        sigi      = idList(idx).freq;
        freqi = obs.range.freqs(obs.range.PRN == prni & obs.range.constInds == constIndi & obs.range.sig == sigi & obs.range.subtype == measTypei);
        
        losInd = losInds(idx);
        weighti = 1 ./ (sin(el(losInd)).^2);
        
        switch measTypei
            case navsu.internal.MeasEnum.Code
                % Code phase measurement
                [predMeasi,Hii,sigMeasi] = codeModel(obj,SimpleModel,nState,sigi,freqi,tecSlant(losInd),x_est_propagated,...
                    constIndi,indGloDcbs(idx),indMpCodes(idx),m(losInd),gRange(losInd),satBias(losInd),rxBias(losInd),trop(losInd),stRangeOffset(losInd),...
                    relClockCorr(losInd),relRangeCorr(losInd),A(losInd,:));                
                
            case navsu.internal.MeasEnum.Carrier
                % Carrier phase measurement
                [predMeasi,Hii,sigMeasi] = carrierModel(obj,nState,sigi,freqi,...
                    tecSlant(losInd),x_est_propagated,m(losInd),indIonos(idx), ...
                    indMpCarrs(idx),indAmbStates(idx),phWind(losInd),...
                    gRange(losInd),satBias(losInd),rxBias(losInd),trop(losInd),...
                    stRangeOffset(losInd),relClockCorr(losInd),relRangeCorr(losInd),...
                    A(losInd,:),constIndi);
                
            case navsu.internal.MeasEnum.Doppler
                % Doppler measurement
                [predMeasi,Hii,sigMeasi] = obj.doppModel(nState,dVel(losInd,:),...
                    A(losInd,:),rxDrift(losInd),constIndi);
                
        end
        
        % Put computed predicted measurement, sensitivity, and
        % measurement sigma into their respective places.
        predMeas(idx) = predMeasi;
        H(idx,:) = Hii;
        R(idx,idx) = weighti*sigMeasi^2;
    end
    
    % Remove measurements from satellites that are too low
    indsElLow = find(el(losInds)*180/pi < PARAMS.elMask);
    if ~isempty(indsElLow)        
        measIdRemovedLow = idList(indsElLow);
        
        H(indsElLow,:) = [];
        R(indsElLow,:) = [];
        R(:,indsElLow) = [];
        predMeas(indsElLow) = [];
        idList(indsElLow) = [];
        measList(indsElLow) = [];
        
    else
        measIdRemovedLow = [];
    end
else
    % No GNSS measurements
    % build each measurement
    predMeas = zeros(0,1);
    H         = zeros(0,nState);
    R         = zeros(0,1);
    measList  = zeros(0,1);
end

if nMeas  == 0
    el = [];
    az = [];
    prnConstInds = zeros(0,2);
    measIdRemovedLow = [];
    idList = [];
end



end