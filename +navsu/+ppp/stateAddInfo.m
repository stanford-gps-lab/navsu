function [measInfoAvail,stateNotNeeded] = stateAddInfo(stateType,gnssMeas,ekf)
% stateAddInfo
% Given the state that is being examined, the currently available  measurements,
% and the current states being tracked, determine what, if any, new states
% need to be added and what can be removed.
% INPUTS:
%   stateType   One of the following:
%                   'cp'           1.  carrier phase
%                   'L1DELAYSTATE' 2. ionospheric LOS delay state
%                   'RX_DCB_GLO'   3. GLONASS receiver differential code bias
%                   'RX_DCB_GPS'   3. GPS receiver diffential code bias
%                   'MP_CODE'      4. code multipath
%                   'MP_CARR'      5. carrier multipath
%   gnssMeas    GNSS observation structure for this epoch- may or may not
%               be needed.
%   ekf         Kalman filter object
% OUPUTS:
%   measInfoAvail  Nx4 vector describing new states to be added. First field
%               is PRN, second field is constellation index, third field is
%               stateType number code (listed above), fourth field is typically
%               the signal number (1,2,3, 102,103 indicating single or dual
%               frequency etc)
%  stateNOtNeeded Just like measInfoAvail except for states that are to be
%               removed.


% Initialize outputs :)
measInfoAvail = [];
stateNotNeeded = [];

switch stateType
    case 'cp'
        % 1
        % Check what carrier phase measurements are currently available
        indsCpAvail = find(gnssMeas.range.ind == 2 & gnssMeas.range.obs > 0);
        
        % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsCpAvail) gnssMeas.range.constInds(indsCpAvail) ...
            ones(size(indsCpAvail)) gnssMeas.range.sig(indsCpAvail)];
        
        stateNotNeeded = find(~ismember(ekf.INDS_STATE.FLEX_STATES_INFO(:,[1 2 4]),measInfoAvail(:,[1 2 4]),'rows') ...
            & ekf.INDS_STATE.FLEX_STATES_INFO(:,3) == 1);
        
        % Look for these in the currently built state
        stateFound = find(ismember(measInfoAvail,ekf.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows'));
        
        % Remove things from the output that were already found in the
        % existing state
        measInfoAvail(stateFound,:) = [];
        
    case 'L1DELAYSTATE'
        % 2
        indsRangeAvail = find( gnssMeas.range.obs > 0);
        
        % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];
        
        prnConstIndsSf = unique(measInfoAvail(measInfoAvail(:,4) < 100,[1 2]),'rows');
        
        % PRN | CONST | STATE TYPE (2) | SIGNAL NUMBER
        measInfoAvail = [prnConstIndsSf 2*ones(size(prnConstIndsSf,1),1) nan(size(prnConstIndsSf,1),1)];
        
        % Remove stuff that is no longer necessary
        stateNotNeeded = find(~ismember(ekf.INDS_STATE.FLEX_STATES_INFO(:,[1 2 3]),...
            [prnConstIndsSf 2*ones(size(prnConstIndsSf,1),1)],'rows') & ekf.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
        
        % Add states that are necessary
        stateFound = ismember(measInfoAvail(:,1:3),ekf.INDS_STATE.FLEX_STATES_INFO(:,1:3),'rows');
        
        measInfoAvail(stateFound,:) = [];
        
    case 'TECSTATE'
        % 2
        indsRangeAvail = find( gnssMeas.range.obs > 0);
        
        % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];
        
        prnConstIndsSf = unique(measInfoAvail(measInfoAvail(:,4) < 100,[1 2]),'rows');
        
        % PRN | CONST | STATE TYPE (2) | SIGNAL NUMBER
        measInfoAvail = [prnConstIndsSf 2*ones(size(prnConstIndsSf,1),1) nan(size(prnConstIndsSf,1),1)];
        
        % Remove stuff that is no longer necessary
        stateNotNeeded = find(~ismember(ekf.INDS_STATE.FLEX_STATES_INFO(:,[1 2 3]),...
            [prnConstIndsSf 2*ones(size(prnConstIndsSf,1),1)],'rows') & ekf.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
        
        % Add states that are necessary
        stateFound = ismember(measInfoAvail(:,1:3),ekf.INDS_STATE.FLEX_STATES_INFO(:,1:3),'rows');
        
        measInfoAvail(stateFound,:) = [];
        
        
    case 'MP_CODE'
        % 4
        % Check what carrier phase measurements are currently available
        indsRangeAvail = find( gnssMeas.range.obs > 0 & gnssMeas.range.ind == 1);
        
        % PRN | CONST | STATE TYPE (4) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) 4*ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];
        
        
        % Remove stuff that is no longer necessary
        stateNotNeeded = find(~ismember(ekf.INDS_STATE.FLEX_STATES_INFO(:,:),...
            measInfoAvail,'rows') & ekf.INDS_STATE.FLEX_STATES_INFO(:,3) == 4);
        
        stateFound = ismember(measInfoAvail,ekf.INDS_STATE.FLEX_STATES_INFO(:,:),'rows');
        
        measInfoAvail(stateFound,:) = [];
        
        
    case 'MP_CARR'
        % 5
        % Check what carrier phase measurements are currently available
        indsRangeAvail = find( gnssMeas.range.obs > 0 & gnssMeas.range.ind == 2);
        
        % PRN | CONST | STATE TYPE (5) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) 5*ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];
        
        % Remove stuff that is no longer necessary
        stateNotNeeded = find(~ismember(ekf.INDS_STATE.FLEX_STATES_INFO(:,:),...
            measInfoAvail,'rows') & ekf.INDS_STATE.FLEX_STATES_INFO(:,3) == 5);
        
        stateFound = ismember(measInfoAvail,ekf.INDS_STATE.FLEX_STATES_INFO(:,:),'rows');
        
        measInfoAvail(stateFound,:) = [];
        
        
    case 'RX_DCB_GLO'
        % 3
        % Check what carrier phase measurements are currently available
        indsGloPr = find(gnssMeas.range.ind == 1 & gnssMeas.range.obs > 0 & ...
            gnssMeas.range.constInds == 2);
        
        % PRN | CONST | STATE TYPE (3) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsGloPr) 2*ones(size(indsGloPr)) ...
            3*ones(size(indsGloPr)) gnssMeas.range.sig(indsGloPr)];
        
        if ~isfield(ekf.INDS_STATE,'RX_DCB_GLO_INFO')
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
            
            ekf.INDS_STATE.RX_DCB_GLO_INFO = [prnHighSnr sigRef];
        end
        % Add states that are necessary
        stateFound = find(ismember(measInfoAvail(:,:),ekf.INDS_STATE.FLEX_STATES_INFO(:,:),'rows') | ...
            (measInfoAvail(:,1) == ekf.INDS_STATE.RX_DCB_GLO_INFO(1) & ...
            measInfoAvail(:,4) == ekf.INDS_STATE.RX_DCB_GLO_INFO(2) ));
        
        measInfoAvail(stateFound,:) = [];
        
    case 'RX_DCB_GPS'
        % 3
        % Check what carrier phase measurements are currently available
        indsGpsPr = find(gnssMeas.range.ind == 1 & gnssMeas.range.obs > 0 & ...
            gnssMeas.range.constInds == 1);
        
        % PRN | CONST | STATE TYPE (3) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsGpsPr) 1*ones(size(indsGpsPr)) ...
            3*ones(size(indsGpsPr)) gnssMeas.range.sig(indsGpsPr)];
        
        % Add states that are necessary
        stateFound = ismember(measInfoAvail(:,:),ekf.INDS_STATE.FLEX_STATES_INFO(:,:),'rows');
        
        if ~isfield(ekf.INDS_STATE,'RX_DCB_GPS_INFO')
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
            
            ekf.INDS_STATE.RX_DCB_GPS_INFO = [prnHighSnr sigRef];
            
        end
        
        stateFound = find(ismember(measInfoAvail(:,:),ekf.INDS_STATE.FLEX_STATES_INFO(:,:),'rows') | ...
            (measInfoAvail(:,1) == ekf.INDS_STATE.RX_DCB_GPS_INFO(1) & ...
            measInfoAvail(:,4) == ekf.INDS_STATE.RX_DCB_GPS_INFO(2) ));
        
        measInfoAvail(stateFound,:) = [];
        
    case 'EPH'
        % 2
        indsRangeAvail = find( gnssMeas.range.obs > 0);
        
        % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsRangeAvail) gnssMeas.range.constInds(indsRangeAvail) ones(size(indsRangeAvail)) gnssMeas.range.sig(indsRangeAvail)];
        
        prnConstInds = unique(measInfoAvail(:,[1 2]),'rows');
        
        % PRN | CONST | STATE TYPE (6) | SIGNAL NUMBER
        measInfoAvail = [prnConstInds 6*ones(size(prnConstInds,1),1) nan(size(prnConstInds,1),1)];
        
        % Remove stuff that is no longer necessary
        stateNotNeeded = find(~ismember(ekf.INDS_STATE.FLEX_STATES_INFO(:,[1 2 3]),...
            [prnConstInds 6*ones(size(prnConstInds,1),1)],'rows') & ekf.INDS_STATE.FLEX_STATES_INFO(:,3) == 6);
        
        % Add states that are necessary
        stateFound = ismember(measInfoAvail(:,1:3),ekf.INDS_STATE.FLEX_STATES_INFO(:,1:3),'rows');
        
        measInfoAvail(stateFound,:) = [];
        
end