function measRemovedSlip = manageStates(obj,epoch,gnssMeas,PARAMS,outStruc)


% Check what carrier phase measurements are currently available
indsCpAvail = find(gnssMeas.range.ind == 2 & gnssMeas.range.obs > 0);

% PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
measInfoAvail = [gnssMeas.range.PRN(indsCpAvail) gnssMeas.range.constInds(indsCpAvail) ones(size(indsCpAvail)) gnssMeas.range.sig(indsCpAvail)];

% State types:
% 1 = carrier phase ambiguity- one per CP
% 2 = TEC state- one per LOS
% 3 = Satellite specific DCB
% 4 = Code phase multipath
% 5 = Carr phase multipath

% Check for cycle slips so that these can be removed and reset
measRemovedSlip = obj.checkCycleSlips(epoch,gnssMeas,PARAMS);

if ~isempty(measRemovedSlip) && ~isempty(outStruc)
  outStruc.saveMeasRemoved(epoch,zeros(0,2),zeros(0,6),measRemovedSlip);
end

% Look for stuff to remove- these are state without associated measurements
stateNotNeeded = find(~ismember(obj.INDS_STATE.FLEX_STATES_INFO(:,[1 2 4]),measInfoAvail(:,[1 2 4]),'rows') ...
    & obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 1);

if ~isempty(stateNotNeeded)
    % Remove from the state estimate, covariance, FLEX_STATES, and
    % FLEX_STATES_INFO
    obj.state(obj.INDS_STATE.FLEX_STATES(stateNotNeeded)) = [];
    obj.cov(obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
    obj.cov(:,obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
    obj.INDS_STATE.FLEX_STATES_INFO(stateNotNeeded,:) = [];
    obj.INDS_STATE.FLEX_STATES = obj.INDS_STATE.FLEX_STATE_MIN-1+...
        [1:size(obj.INDS_STATE.FLEX_STATES_INFO,1)]';
end

% Look for these in the currently built state
stateFound = ismember(measInfoAvail,obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');

for idx = 1:length(indsCpAvail)
    if stateFound(idx)
        continue;
    end
    
    % Add this new state
    infoAdd = measInfoAvail(idx,:);
    
    switch infoAdd(3)
        case 1
            % Carrier phase
            
            % Produce an initial carrier phase estimate
            % Find the corresponding code phase estimate and just match to
            % that
            indCode = find(gnssMeas.range.ind == 1 & gnssMeas.range.PRN == infoAdd(1) & ...
                gnssMeas.range.constInds == infoAdd(2) & gnssMeas.range.sig == infoAdd(4));
            indCp = find(gnssMeas.range.ind == 2 & gnssMeas.range.PRN == infoAdd(1) & ...
                gnssMeas.range.constInds == infoAdd(2) & gnssMeas.range.sig == infoAdd(4));
            
            obsCodePhase = gnssMeas.range.obs(indCode);
            obsCarrPhase = gnssMeas.range.obs(indCp);
            
            if obsCodePhase == 0
                ambEst0 = 0;
            else
                ambEst0 = -obsCodePhase+obsCarrPhase;
            end
            
            ambCov0 = PARAMS.SIGMA0.AMB;
            
            indStateAdd = obj.INDS_STATE.FLEX_STATE_MIN+length(obj.INDS_STATE.FLEX_STATES);
            
            % Add to the state, covariance, FLEX_STATE_INFO, and
            % FLEX_STATES objects
            
            % state
            obj.state = [obj.state; ambEst0];
            
            % covariance
            covNew = zeros(indStateAdd);
            covNew(1:size(obj.cov,1),1:size(obj.cov,1)) = obj.cov;
            covNew(indStateAdd,indStateAdd) = ambCov0^2;
            obj.cov = covNew;
            
            % FLEX_STATE_INFO
            obj.INDS_STATE.FLEX_STATES_INFO = [obj.INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];

            % FLEX_STATE index
            obj.INDS_STATE.FLEX_STATES = [obj.INDS_STATE.FLEX_STATES; indStateAdd];
    end
end

%% TEC state
if PARAMS.states.iono && strcmp(PARAMS.states.ionoMode,'L1DELAYSTATE')
    % Check what carrier phase measurements are currently available
    indsRangeAvail = find( gnssMeas.range.obs > 0);
    
    % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
    measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];

    prnConstIndsSf = unique(measInfoAvail(measInfoAvail(:,4) < 100,[1 2]),'rows');
    
    stateInfoIono = [prnConstIndsSf 2*ones(size(prnConstIndsSf,1),1) nan(size(prnConstIndsSf,1),1)];
    
    % Remove stuff that is no longer necessary
    stateNotNeeded = find(~ismember(obj.INDS_STATE.FLEX_STATES_INFO(:,[1 2 3]),...
        [prnConstIndsSf 2*ones(size(prnConstIndsSf,1),1)],'rows') & obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
    if ~isempty(stateNotNeeded)
        % Remove from the state estimate, covariance, FLEX_STATES, and
        % FLEX_STATES_INFO
        obj.state(obj.INDS_STATE.FLEX_STATES(stateNotNeeded)) = [];
        obj.cov(obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.cov(:,obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.INDS_STATE.FLEX_STATES_INFO(stateNotNeeded,:) = [];
        obj.INDS_STATE.FLEX_STATES = obj.INDS_STATE.FLEX_STATE_MIN-1+...
            [1:size(obj.INDS_STATE.FLEX_STATES_INFO,1)]';
    end
    
    % Add states that are necessary
    stateFound = ismember(stateInfoIono(:,1:3),obj.INDS_STATE.FLEX_STATES_INFO(:,1:3),'rows');
    
    for idx = 1:length(stateFound)
        if stateFound(idx)
            continue;
        end
        
        infoAdd = stateInfoIono(idx,:);
        
        indStateAdd = obj.INDS_STATE.FLEX_STATE_MIN+length(obj.INDS_STATE.FLEX_STATES);
        
        % Add to the state, covariance, FLEX_STATE_INFO, and
        % FLEX_STATES objects
        
        % state
        obj.state = [obj.state; 0];
        
        % covariance
        covNew = zeros(indStateAdd);
        covNew(1:size(obj.cov,1),1:size(obj.cov,1)) = obj.cov;
        covNew(indStateAdd,indStateAdd) = PARAMS.SIGMA0.L1_IONO^2;
        obj.cov = covNew;
        
        % FLEX_STATE_INFO
        obj.INDS_STATE.FLEX_STATES_INFO = [obj.INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];
        
        % FLEX_STATE index
        obj.INDS_STATE.FLEX_STATES = [obj.INDS_STATE.FLEX_STATES; indStateAdd];
    end    
end

%% GLONASS Satellite-Receiver DCB
% One DCB per satellite and signal- accounts for differing delays per
% frequency in the receiver
if PARAMS.states.RX_DCB_GLO
    % Check what carrier phase measurements are currently available
    indsGloPr = find(gnssMeas.range.ind == 1 & gnssMeas.range.obs > 0 & ...
        gnssMeas.range.constInds == 2);
    
    % GLONASS PRNs and signals 
    stateInfoGloIono = [gnssMeas.range.PRN(indsGloPr) 2*ones(size(indsGloPr)) ...
        3*ones(size(indsGloPr)) gnssMeas.range.sig(indsGloPr)];
    
    % Add states that are necessary
    stateFound = ismember(stateInfoGloIono(:,:),obj.INDS_STATE.FLEX_STATES_INFO(:,:),'rows');
    
    if ~isfield(obj.INDS_STATE,'RX_DCB_GLO_INFO')
        % choose the signal with the highest SNR
        
        indsGlo = find(gnssMeas.snr.constInds == 2 & gnssMeas.snr.sig == 1);
        [~,indHighSnr] = max(gnssMeas.snr.obs(indsGlo));
        prnHighSnr = gnssMeas.snr.PRN(indsGlo(indHighSnr));
        snrHigh = gnssMeas.snr.obs(indsGlo(indHighSnr));
        
        indsPrnHigh = find(gnssMeas.snr.PRN == prnHighSnr & gnssMeas.snr.constInds == 2 & gnssMeas.snr.sig ~= 1);
        
        [~,indHighSnr2] = max(gnssMeas.snr.obs(indsPrnHigh));
        snrHigh2 = gnssMeas.snr.obs(indsPrnHigh(indHighSnr2));
        
        if snrHigh2 > 0
            % signal is dual frequency combination
            sig2 = gnssMeas.snr.sig(indsPrnHigh(indHighSnr2));
            
            sigRef = 100+sig2;
        else
            sigRef = 1;
        end
        
        obj.INDS_STATE.RX_DCB_GLO_INFO = [prnHighSnr sigRef];
    end
    
    for idx = 1:length(stateFound)
        
        if stateFound(idx) || (stateInfoGloIono(idx,1) == obj.INDS_STATE.RX_DCB_GLO_INFO(1) && ...
                stateInfoGloIono(idx,4) == obj.INDS_STATE.RX_DCB_GLO_INFO(2) )
            continue;
        end
        
        %         if stateFound(idx) || (stateInfoGloIono(idx,1) == obj.INDS_STATE.RX_DCB_GLO_INFO(1))
        %             continue;
        %         end
        
        infoAdd = stateInfoGloIono(idx,:);
        
        indStateAdd = obj.INDS_STATE.FLEX_STATE_MIN+length(obj.INDS_STATE.FLEX_STATES);
        
        % Add to the state, covariance, FLEX_STATE_INFO, and
        % FLEX_STATES objects
        
        % state
        obj.state = [obj.state; 0];
        
        % covariance
        covNew = zeros(indStateAdd);
        covNew(1:size(obj.cov,1),1:size(obj.cov,1)) = obj.cov;
        covNew(indStateAdd,indStateAdd) = PARAMS.SIGMA0.RX_DCB_GLO.^2;
        obj.cov = covNew;
        
        % FLEX_STATE_INFO
        obj.INDS_STATE.FLEX_STATES_INFO = [obj.INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];
        
        % FLEX_STATE index
        obj.INDS_STATE.FLEX_STATES = [obj.INDS_STATE.FLEX_STATES; indStateAdd];
    end
end


%% GPS Satellite-Receiver DCB
% One DCB per satellite and signal- accounts for differing delays per
% frequency in the receiver

if PARAMS.states.RX_DCB_GPS
    % Check what carrier phase measurements are currently available
    indsGpsPr = find(gnssMeas.range.ind == 1 & gnssMeas.range.obs > 0 & ...
        gnssMeas.range.constInds == 1);
    
    % GLONASS PRNs and signals 
    stateInfoGpsIono = [gnssMeas.range.PRN(indsGpsPr) 1*ones(size(indsGpsPr)) ...
        3*ones(size(indsGpsPr)) gnssMeas.range.sig(indsGpsPr)];
    
    % Add states that are necessary
    stateFound = ismember(stateInfoGpsIono(:,:),obj.INDS_STATE.FLEX_STATES_INFO(:,:),'rows');
    
    if ~isfield(obj.INDS_STATE,'RX_DCB_GPS_INFO')
        % choose the signal with the highest SNR
        
        indsGlo = find(gnssMeas.snr.constInds == 1 & gnssMeas.snr.sig == 1);
        [~,indHighSnr] = max(gnssMeas.snr.obs(indsGlo));
        prnHighSnr = gnssMeas.snr.PRN(indsGlo(indHighSnr));
        snrHigh = gnssMeas.snr.obs(indsGlo(indHighSnr));
        
        indsPrnHigh = find(gnssMeas.snr.PRN == prnHighSnr & gnssMeas.snr.constInds == 1 & gnssMeas.snr.sig ~= 1);
        
        [~,indHighSnr2] = max(gnssMeas.snr.obs(indsPrnHigh));
        snrHigh2 = gnssMeas.snr.obs(indsPrnHigh(indHighSnr2));
        
        if snrHigh2 > 0
            % signal is dual frequency combination
            sig2 = gnssMeas.snr.sig(indsPrnHigh(indHighSnr2));
            
            sigRef = 100+sig2;
        else
            sigRef = 1;
        end
        
        obj.INDS_STATE.RX_DCB_GPS_INFO = [prnHighSnr sigRef];
    end
    
    for idx = 1:length(stateFound)
        
        if stateFound(idx) || (stateInfoGpsIono(idx,1) == obj.INDS_STATE.RX_DCB_GPS_INFO(1) && ...
                stateInfoGpsIono(idx,4) == obj.INDS_STATE.RX_DCB_GPS_INFO(2) )
            continue;
        end
        
        %         if stateFound(idx) || (stateInfoGpsIono(idx,1) == obj.INDS_STATE.RX_DCB_GLO_INFO(1))
        %             continue;
        %         end
        
        infoAdd = stateInfoGpsIono(idx,:);
        
        indStateAdd = obj.INDS_STATE.FLEX_STATE_MIN+length(obj.INDS_STATE.FLEX_STATES);
        
        % Add to the state, covariance, FLEX_STATE_INFO, and
        % FLEX_STATES objects
        
        % state
        obj.state = [obj.state; 0];
        
        % covariance
        covNew = zeros(indStateAdd);
        covNew(1:size(obj.cov,1),1:size(obj.cov,1)) = obj.cov;
        covNew(indStateAdd,indStateAdd) = PARAMS.SIGMA0.RX_DCB_GPS.^2;
        obj.cov = covNew;
        
        % FLEX_STATE_INFO
        obj.INDS_STATE.FLEX_STATES_INFO = [obj.INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];
        
        % FLEX_STATE index
        obj.INDS_STATE.FLEX_STATES = [obj.INDS_STATE.FLEX_STATES; indStateAdd];
    end
    
end


%% Code phase multipath
if PARAMS.states.MP_CODE
    % Check what carrier phase measurements are currently available
    indsRangeAvail = find( gnssMeas.range.obs > 0 & gnssMeas.range.ind == 1);
    
    % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
    measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) 4*ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];

%     prnConstIndsSf = unique(measInfoAvail(measInfoAvail(:,4) < 100,[1 2]),'rows');
    
    stateInfoIono = measInfoAvail;
    
    % Remove stuff that is no longer necessary
    stateNotNeeded = find(~ismember(obj.INDS_STATE.FLEX_STATES_INFO(:,:),...
        measInfoAvail,'rows') & obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 4);
    if ~isempty(stateNotNeeded)
        % Remove from the state estimate, covariance, FLEX_STATES, and
        % FLEX_STATES_INFO
        obj.state(obj.INDS_STATE.FLEX_STATES(stateNotNeeded)) = [];
        obj.cov(obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.cov(:,obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.INDS_STATE.FLEX_STATES_INFO(stateNotNeeded,:) = [];
        obj.INDS_STATE.FLEX_STATES = obj.INDS_STATE.FLEX_STATE_MIN-1+...
            [1:size(obj.INDS_STATE.FLEX_STATES_INFO,1)]'; 
    end
    
    % Add states that are necessary
    stateFound = ismember(measInfoAvail,obj.INDS_STATE.FLEX_STATES_INFO(:,:),'rows');
    
    for idx = 1:length(stateFound)
        if stateFound(idx)
            continue;
        end
        
        infoAdd = stateInfoIono(idx,:);
        
        indStateAdd = obj.INDS_STATE.FLEX_STATE_MIN+length(obj.INDS_STATE.FLEX_STATES);
        
        % Add to the state, covariance, FLEX_STATE_INFO, and
        % FLEX_STATES objects
        
        % state
        obj.state = [obj.state; 0];
        
        % covariance
        covNew = zeros(indStateAdd);
        covNew(1:size(obj.cov,1),1:size(obj.cov,1)) = obj.cov;
        covNew(indStateAdd,indStateAdd) = PARAMS.SIGMA0.MP_CODE^2;
        obj.cov = covNew;
        
        % FLEX_STATE_INFO
        obj.INDS_STATE.FLEX_STATES_INFO = [obj.INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];
        
        % FLEX_STATE index
        obj.INDS_STATE.FLEX_STATES = [obj.INDS_STATE.FLEX_STATES; indStateAdd];
    end
end

%% Carrier phase multipath
if PARAMS.states.MP_CARR
    % Check what carrier phase measurements are currently available
    indsRangeAvail = find( gnssMeas.range.obs > 0 & gnssMeas.range.ind == 2);
    
    % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
    measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) 4*ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];

%     prnConstIndsSf = unique(measInfoAvail(measInfoAvail(:,4) < 100,[1 2]),'rows');
    
    stateInfoIono = measInfoAvail;
    
    % Remove stuff that is no longer necessary
    stateNotNeeded = find(~ismember(obj.INDS_STATE.FLEX_STATES_INFO(:,:),...
        measInfoAvail,'rows') & obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 5);
    if ~isempty(stateNotNeeded)
        % Remove from the state estimate, covariance, FLEX_STATES, and
        % FLEX_STATES_INFO
        obj.state(obj.INDS_STATE.FLEX_STATES(stateNotNeeded)) = [];
        obj.cov(obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.cov(:,obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.INDS_STATE.FLEX_STATES_INFO(stateNotNeeded,:) = [];
        obj.INDS_STATE.FLEX_STATES = obj.INDS_STATE.FLEX_STATE_MIN-1+...
            [1:size(obj.INDS_STATE.FLEX_STATES_INFO,1)]'; 
    end
    
    % Add states that are necessary
    stateFound = ismember(measInfoAvail,obj.INDS_STATE.FLEX_STATES_INFO(:,:),'rows');
    
    for idx = 1:length(stateFound)
        if stateFound(idx)
            continue;
        end
        
        infoAdd = stateInfoIono(idx,:);
        
        indStateAdd = obj.INDS_STATE.FLEX_STATE_MIN+length(obj.INDS_STATE.FLEX_STATES);
        
        % Add to the state, covariance, FLEX_STATE_INFO, and
        % FLEX_STATES objects
        
        % state
        obj.state = [obj.state; 0];
        
        % covariance
        covNew = zeros(indStateAdd);
        covNew(1:size(obj.cov,1),1:size(obj.cov,1)) = obj.cov;
        covNew(indStateAdd,indStateAdd) = PARAMS.SIGMA0.MP_CODE^2;
        obj.cov = covNew;
        
        % FLEX_STATE_INFO
        obj.INDS_STATE.FLEX_STATES_INFO = [obj.INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];
        
        % FLEX_STATE index
        obj.INDS_STATE.FLEX_STATES = [obj.INDS_STATE.FLEX_STATES; indStateAdd];
    end
end

end