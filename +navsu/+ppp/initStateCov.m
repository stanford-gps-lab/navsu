function [state0,cov0] = initStateCov(stateType,infoAdd,PARAMS,gnssMeas)
% initStateCov
% Produces an initial estimate and covariance for various "flex states," 
% which can be any several things. 
% INPUTS:
%   stateType   One of the following: 
%                   'cp'           1.  carrier phase
%                   'L1DELAYSTATE' 2. ionospheric LOS delay state
%                   'RX_DCB_GLO'   3. GLONASS receiver differential code bias
%                   'RX_DCB_GPS'   3. GPS receiver diffential code bias
%                   'MP_CODE'      4. code multipath
%                   'MP_CARR'      5. carrier multipath
%                   'EPH'          6. ephemeris 
%   infoAdd     1x4 vector. First field is PRN, second field is
%               constellation index, third field is stateType number code
%               (listed above), fourth field is typically the signal number
%               (1,2,3, 102,103 indicating single or dual frequency etc)
%   PARAMS      typical PPP PARAMS object
%   gnssMeas    GNSS observation structure for this epoch- may or may not
%               be needed.
% OUTPUTS:
%   state0      initial state estimate
%   cov0        initial covariance for this estimate
state0 = [];
cov0 = [];

switch stateType
    case 'cp'
        % Carrier phase ambiguity estimate
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
        
        % state
        state0 = ambEst0;
        
        % covariance
        cov0 = PARAMS.SIGMA0.AMB^2;
        
    case 'L1DELAYSTATE'
        % Ionospheric delay state- delay on L1 per LOS
        state0 = 0;
        
        cov0 =  PARAMS.SIGMA0.L1_IONO^2;
        
    case 'TECSTATE'
        % Ionospheric delay state that uses a TEC map for initial
        % correction on single frequency measurements
        state0 = 0;
        
        cov0 =  PARAMS.SIGMA0.L1_IONO^2;
        
    case 'MP_CODE'
        % Code multipath
        state0 = 0;
        cov0 = PARAMS.SIGMA0.MP_CODE^2;
        
    case 'MP_CARR'
        % Carrier multipath
        state0 = 0;
        cov0   =  PARAMS.SIGMA0.MP_CARR^2;
        
        
    case 'RX_DCB_GLO'
        state0 = 0;
        cov0   = PARAMS.SIGMA0.RX_DCB_GLO.^2;
        
    case 'RX_DCB_GPS'
        state0 = 0;
        cov0   = PARAMS.SIGMA0.RX_DCB_GPS.^2;
        
    case 'EPH' 
        state0 = 0;
        cov0 = 2.4^2;
end


















end