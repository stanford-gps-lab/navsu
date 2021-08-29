function [state,dstate,resids,measMat,covState,covdState,satsUsed] = lsSolGnss(obs,corrData,PARAMS,varargin)

%% WARNING function not used throughout PPP filters, could be faulty.

% Parse optional inputs
p = inputParser;
p.addParameter('pos0', [0 0 0]); % whe
parse(p, varargin{:});
res = p.Results;
pos0 = res.pos0;

% Unique constellations
constUn = unique(obs.constInds);
nConst0  = length(constUn);

covState  = [];
covdState = [];
satsUsed  = [];
%% Collect all of the measurements we want- prioritize dual frequency
% indsDfCodeAvail = find(obs.range.obs ~= 0 & obs.range.sig >= 100 & ...
%     obs.range.ind == 1);

indsDfCodeAvail = find(obs.range.obs(:) ~= 0 & [obs.range.ID.freq]' >= 100 & ...
    [obs.range.ID.subtype]' == navsu.internal.MeasEnum.Code);

% Also include single frequency if we have to :)
if size(indsDfCodeAvail,1) < 1000
    indsSfCodeAvail = find(obs.range.obs(:) ~= 0 &  [obs.range.ID.freq]' < 100 & ...
        [obs.range.ID.subtype]' ==  navsu.internal.MeasEnum.Code);
else
    indsSfCodeAvail = [];
end

indsMeasUse = [indsSfCodeAvail; indsDfCodeAvail];

prnObsMat      = repmat(obs.PRN,size(obs.range.obs,1),1);
constIndObsMat = repmat(obs.constInds,size(obs.range.obs,1),1);
% freqDopMat     = repmat(obs.freqs(1,:)

nMeas      = length(indsMeasUse);
% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas
measMat = [prnObsMat(indsMeasUse) constIndObsMat(indsMeasUse) ...
    obs.range.sig(indsMeasUse) obs.freqs(indsMeasUse) obs.range.obs(indsMeasUse) ];
pr = measMat(:,5);

%% Collect doppler measurements
indsMeasUseDop = find(obs.doppler.obs ~= 0) ;

prnDopMat      = repmat(obs.PRN,size(obs.doppler.obs,1),1);
constIndDopMat = repmat(obs.constInds,size(obs.doppler.obs,1),1);

nMeasDop      = length(indsMeasUseDop);
% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas
measMatDop = [prnDopMat(indsMeasUseDop) constIndDopMat(indsMeasUseDop) ...
    obs.doppler.sig(indsMeasUseDop) nan(size(indsMeasUseDop)) obs.doppler.obs(indsMeasUseDop) ];

%% 
epoch = obs.epochs(1);

c = navsu.constants.c;

% Unique constellations
nConst  = length(constUn);

% Number of least squares iterations...
nIter = 10;

% Initialize position and clock state
state = zeros(3+nConst,1);
state(1:3) = pos0;
dstate = zeros(3+nConst,1);

rxBias = zeros(nMeas,1);
trop0  = zeros(nMeas,1);

% Satellite clock bias
sbias = c*corrData.clock(measMat(:,1),measMat(:,2),epoch*ones(nMeas,1));

% Initialize iono corrections
ionoDelay = zeros(nMeas,1);

% Receiver clock bias
rbias = zeros(nMeas,1);

% Initialize signal transmission time
ttx = epoch-pr/c+rbias+trop0/c-sbias/c*0-ionoDelay;

% Initialize satellite positions
[svPos,svVel] = corrData.propagate(measMat(:,1),measMat(:,2),obs.epochs(1)*ones(nMeas,1));
    
% If there are any satellites with missing precise data, remove them here
rdx = find(isnan(svPos(:,1)));
svPos(rdx,:)     = [];
svVel(rdx,:)     = [];
pr(rdx)          = [];
measMat(rdx,:)   = [];
trop0(rdx)       = [];
sbias(rdx)       = [];
indsMeasUse(rdx) = [];
rxBias(rdx)      = [];
ionoDelay(rdx)   = [];
% doppler(rdx)    = [];

nMeas = nMeas-length(rdx);

% Compute day of year for tropo computation
doy = navsu.time.jd2doy(navsu.time.epochs2jd(epoch));
rxBiasByObs = rxBias;
epochi = epoch;
xNorm = 1000;
closeThresh = 500;
for idx = 1:nIter
    usrPos = state(1:3);
    usrVel = dstate(1:3);
    
    epochi = epoch-rxBiasByObs/c;

    % Receiver constellation clock bias
    [~,bdx] = ismember(measMat(:,2)',constUn);
    rxBiasByObs = state(bdx+3);
    rxBiasRateByObs = dstate(bdx+3);
    
    ttx = epochi-pr/c+rxBiasByObs/c*0+trop0/c*0+sbias/c*0;
    
    % Remove anything with a NaN
    if any(isnan(ttx))
        pr(isnan(ttx))           = [];
        measMat(isnan(ttx),:)    = [];
        sbias(isnan(ttx))        = [];
        svPos(isnan(ttx),:)      = [];
        svVel(isnan(ttx),:)      = [];
        indsMeasUse(isnan(ttx))  = [];        
        ttx(isnan(ttx))          = [];
        ionoDelay(isnan(ttx))    = [];
        
        nMeas = length(pr);
        [~,bdx] = ismember(measMat(:,2)',constUn);
        rxBiasByObs = state(bdx+3);
    end
    
    % Update satellite position
    if xNorm < closeThresh
        [svPos,svVel,~,sigOrbit] = corrData.propagate(measMat(:,1),measMat(:,2),ttx);
    end
    
    travelTime = epochi-ttx;
    svPosRot = navsu.ppp.models.earthRotTravelCorr(travelTime,svPos);
    
    % Tropospheric delay
    if xNorm < closeThresh
        % Compute elevation and azimuth angle
        [el,az] = navsu.geo.pos2elaz(usrPos',svPos);
        
        elMask = -10;
        
        % remove doppler measurements if el is too low
        prnConstRemove = measMat(el < elMask*pi/180,1:2);
        measMatDop(ismember(measMatDop(:,1:2),prnConstRemove,'rows'),:) = [];
        
        pr(el < elMask*pi/180)           = [];
        measMat(el < elMask*pi/180,:)    = [];
        sbias(el < elMask*pi/180)        = [];
        svPos(el < elMask*pi/180,:)      = [];
        sigOrbit(el < elMask*pi/180,:)   = [];
        svVel(el < elMask*pi/180,:)      = [];
        indsMeasUse(el < elMask*pi/180)  = [];
        svPosRot(el < elMask*pi/180,:)   = [];
        ionoDelay(el < elMask*pi/180)    = [];
        ttx(isnan(ttx))                  = [];        
        
        [~,bdx] = ismember(measMat(:,2)',constUn);
        rxBiasByObs = state(bdx+3);
        rxBiasRateByObs = dstate(bdx+3);
        
        az(el < elMask*pi/180) = [];
        el(el < elMask*pi/180) = [];
        
        nMeas = length(pr);
        
        usrLlh = navsu.geo.xyz2llh(usrPos');
        
        % Update tropo delay
        trop0 = navsu.ppp.models.saastamoinen(usrLlh(1), usrLlh(2), usrLlh(3), el*180/pi);
        
        % Update iono delay
        ionoDelay = zeros(nMeas,1);
        indsSf = find(measMat(:,3) < 100);
        if ~isempty(indsSf)
            [~,ionoDelay(indsSf)] = corrData.ionoDelay(epochi(1),usrLlh,'az',az(indsSf),...
                'el',el(indsSf),'satPos',svPos(indsSf,:),'freqs',measMat(indsSf,4));
        end
        
    else
        trop0 = zeros(nMeas,1);
        ionoDelay = zeros(nMeas,1);
    end
    
    % Geometric range
    dPos = svPosRot-repmat(usrPos',nMeas,1);
    gRange = sqrt(sum(dPos.^2,2));
    
    % Relativistic clock correction
    relCorr = 2/c^2.*sum(svPos.*svVel,2)*c;
    prEst = gRange+rxBiasByObs+trop0-sbias+relCorr+ionoDelay;
    resids = pr-prEst;
    
    A = zeros(nMeas,3+nConst);
%     A(:,1:3) = -(svPosRot-repmat(usrPos',nSat,1))./repmat(pr,1,3);
    A(:,1:3) = (repmat(usrPos',nMeas,1)-svPosRot)./sqrt(sum((repmat(usrPos',nMeas,1)-svPosRot).^2,2));
    
    A(sub2ind(size(A),1:size(A,1),bdx+3)) = 1;
    
    if xNorm < closeThresh       
        % Baseline measurement accuracy (un-corrected single freq) 
        Rdiag0 = (PARAMS.sigMeas.iono.^2+PARAMS.sigMeas.pr.^2) *ones(size(el))+sigOrbit.^2;
        % Dual frequency measurement accuracy (conservative!)
        Rdiag0(measMat(:,3) >= 100) = PARAMS.sigMeas.pr.^2+sigOrbit(measMat(:,3) >= 100).^2;
        % Elevation scaling
        elScale = 1./(sin(el).^2);
        W = diag(1./(Rdiag0.*elScale));        
    else
        W = eye(nMeas);
    end
    
    % if we don't have enough measurements, there's no point in doing this
    if length(resids) < 3+nConst
        break;
    end
    
    x =(A'*W*A)\A'*W* resids;
    state = state+x;
    covState = inv(A'*W*A);
    
    % do velocity update
    if xNorm < closeThresh 
        usrPos = state(1:3);
        usrVel = dstate(1:3);
        
        % Sort out what satellite positiions and velocities are available
        svPosRotDop = nan(size(measMatDop,1),3);
        svVelDop    = nan(size(measMatDop,1),3);
        elDop       = nan(size(measMatDop,1),1);
%         [ia,ib] = ismember(measMat(:,1:2),measMatDop(:,1:2),'rows');
        [ia,ib] = ismember(measMatDop(:,1:2),measMat(:,1:2),'rows');
        svPosRotDop(ia,:) = svPosRot(ib(ia),:);
        svVelDop(ia,:)    = svVel(ib(ia),:);
        elDop(ia)         = el(ib(ia));
        
        % remove any entries without positions
        measMatDop(isnan(svPosRotDop(:,1)),:) = [];        
        svVelDop(isnan(svPosRotDop(:,1)),:) = [];
        elDop(isnan(svPosRotDop(:,1)))   = [];
        svPosRotDop(isnan(svPosRotDop(:,1)),:) = [];
        
        
        % Match constellations to measurements
        [~,bdx] = ismember(measMatDop(:,2)',constUn);
        
        nMeasDop = size(measMatDop,1);
        % Build sensitivity matrix
        AVel = zeros(size(measMatDop,1),3+nConst);
        AVel(:,1:3) = -(repmat(usrPos',nMeasDop,1)-svPosRotDop)./sqrt(sum((repmat(usrPos',nMeasDop,1)-svPosRotDop).^2,2));
        AVel(sub2ind(size(AVel),1:size(AVel,1),bdx+3)) = -1;
        
        % Build doppler prediction
        dVel = repmat(usrVel',nMeasDop,1)-svVelDop;
        rxBiasRateByObs = dstate(bdx+3);
        dopplerEst =   dot(dVel',AVel(:,1:3)')' - rxBiasRateByObs;
        
        residsVel = measMatDop(:,5) - dopplerEst;

        % Baseline measurement accuracy (un-corrected single freq) 
        Rdiag0 = PARAMS.sigMeas.dopp.^2*ones(size(elDop))+PARAMS.sigMeas.ionoRate.^2;
        % Dual frequency measurement accuracy 
        Rdiag0(measMatDop(:,3) >= 100) = PARAMS.sigMeas.dopp.^2;
        % Elevation scaling
        elScale = 1./(sin(elDop).^2);
        W = diag(1./(Rdiag0.*elScale));        
        
        xVel =(AVel'*W*AVel)\AVel'*W* residsVel;
        dstate = dstate+xVel;
        covdState = inv(AVel'*W*AVel);
    end    
    
    if xNorm < 1e-1 && norm(xVel) < 5e-1
        break
    end
    xNorm = norm(x);

end

% Compute elevation angles to each satellite
[el,az] = navsu.geo.pos2elaz(usrPos',svPos);

if length(state) < (3+nConst0) 
    % Need to populate all clock states... even if not available
    
    state2 = nan(3+nConst0,1);
    state2(1:3) = state(1:3);
    
    state2(3+constUn) = state(4:end);
    
    state = state2;
end

if nMeas < 4
    state = zeros(size(state));
end
  
% List of satellites used in this solution- 
% prn | constInd
satsUsed = unique(measMat(:,1:2),'rows');

end


























