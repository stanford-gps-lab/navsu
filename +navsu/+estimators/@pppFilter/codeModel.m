function [predMeasi,Hii,sig] = codeModel(obj,SimpleModel,nState,sigi,freqi,tecSlant,state,...
    constIndi,indGloDcbsi,indMpCodesi,m,gRange,satBias,rxBias,trop,stRangeOffset,...
    relClockCorr,relRangeCorr,A,indIonosi,indEphErri)


predMeas = [];
H = [];
sig = [];
Hii = zeros(1,nState);

%% code
if sigi < 100 && strcmp(obj.PARAMS.states.ionoMode,'TEC')
    % If it's single frequency, need to include iono correction
    ionoCorri = -tecSlant*40.3*10^15./freqi.^2;
    
elseif ~SimpleModel && sigi < 100 && strcmp(obj.PARAMS.states.ionoMode,'L1DELAYSTATE')
    
    ionoCorrModel = -tecSlant*40.3*10^15./freqi.^2;
    
    indIono = obj.INDS_STATE.FLEX_STATES(indIonosi);
    delayL1i = state(indIono);
    
    hi = -(1575.42e6).^2./freqi.^2;
    Hii(1,indIono) = hi;
    ionoCorri = delayL1i*hi + ionoCorrModel;
else
    ionoCorri = 0;
end

if ~SimpleModel && obj.PARAMS.states.RX_DCB && ~(obj.PARAMS.states.RX_DCB_GLO && constIndi == 2) ...
        && ~(obj.PARAMS.states.RX_DCB_GPS && constIndi == 1)
    indRxDcb = find(obj.INDS_STATE.RX_DCB.sig == sigi & obj.INDS_STATE.RX_DCB.constInds == constIndi, 1);
    if ~isempty(indRxDcb)
        rxDcb = state(obj.INDS_STATE.RX_DCB.INDS(indRxDcb));
        Hii(1,obj.INDS_STATE.RX_DCB.INDS(indRxDcb)) = 1;
    else
        rxDcb = 0;
    end
else
    rxDcb = 0;
end

if ~SimpleModel && obj.PARAMS.states.RX_DCB_GPS && constIndi == 1 && indGloDcbs(idx)~= 0
    indDcbGlo = obj.INDS_STATE.FLEX_STATES(indGloDcbsi);
    dcbGpsi = x_est_propagated(indDcbGlo);
    Hii(1,indDcbGlo) = 1;
else
    dcbGpsi = 0;
end

if ~SimpleModel && obj.PARAMS.states.RX_DCB_GLO && constIndi == 2 && indGloDcbs(idx)~= 0
    indDcbGlo = obj.INDS_STATE.FLEX_STATES(indGloDcbsi);
    dcbGloi = x_est_propagated(indDcbGlo);
    Hii(1,indDcbGlo) = 1;
else
    dcbGloi = 0;
end

if ~SimpleModel && obj.PARAMS.states.MP_CODE
    indMpCode = obj.INDS_STATE.FLEX_STATES(indMpCodesi);
    mpCodei   = state(indMpCode);
    Hii(1,indMpCode) = 1;
else
    mpCodei = 0;
end

if ~SimpleModel && obj.PARAMS.states.trop
    dtrop = m*state(obj.INDS_STATE.TROP);
    Hii(1,obj.INDS_STATE.TROP) = m;
else
    dtrop = 0;
end

if ~SimpleModel && obj.PARAMS.states.EPH
    indEphErr = obj.INDS_STATE.FLEX_STATES(indEphErri);
    ephErri   = state(indEphErr);
    Hii(1,indEphErr) = 1;
else
    ephErri = 0;
end

predMeasi = gRange+satBias+rxBias+trop+dtrop+...
    stRangeOffset+relClockCorr+relRangeCorr+ionoCorri+rxDcb+...
    dcbGloi+dcbGpsi+ephErri;

Hii(1,obj.INDS_STATE.POS)        = A;
Hii(1,obj.INDS_STATE.CLOCK_BIAS(obj.INDS_STATE.CLOCK_BIAS_CONSTS == constIndi)) = 1;


% Also include lever arm sensitivity lol
% just constructing things

% z axis (
rArmBody = obj.PARAMS.IMU_ARM;
if ~SimpleModel && norm(rArmBody) > 0  
    
    rArmBodyNorm = rArmBody./norm(rArmBody);
    
    % Convert the LOS to the body frame
    rhoBody = obj.R_b_e'*A';
    
    Hangle(1) = sum(norm(rArmBody([2 3])).*cross([1 0 0]',rArmBodyNorm).*rhoBody);
    Hangle(2) = sum(norm(rArmBody([1 3])).*cross([0 1 0]',rArmBodyNorm).*rhoBody);
    Hangle(3) = sum(norm(rArmBody([1 2])).*cross([0 0 1]',rArmBodyNorm).*rhoBody);
    Hii(1,obj.INDS_STATE.ATTITUDE) = Hangle;
    
end
sig = 3;









end