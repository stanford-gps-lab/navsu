function [predMeasi,Hii,sig] = codeModel(obj,nState,sigi,freqi,tecSlant,state,...
    constIndi,indGloDcbsi,indMpCodesi,m,gRange,satBias,rxBias,trop,stRangeOffset,...
    relClockCorr,relRangeCorr,A)


predMeas = [];
H = [];
sig = [];
Hii = zeros(1,nState);

%% code
if sigi < 100 && strcmp(obj.PARAMS.states.ionoMode,'TEC')
    % If it's single frequency, need to include iono correction
    ionoCorri = -tecSlant*40.3*10^15./freqi.^2;
    
elseif sigi < 100 && strcmp(obj.PARAMS.states.ionoMode,'L1DELAYSTATE')
    
    ionoCorrModel = -tecSlant*40.3*10^15./freqi.^2;
    
    indIono = obj.INDS_STATE.FLEX_STATES(indIonos(idx));
    delayL1i = state(indIono);
    
    hi = -(1575.42e6).^2./freqi.^2;
    Hii(1,indIono) = hi;
    ionoCorri = delayL1i*hi+ionoCorrModel;
else
    ionoCorri = 0;
end

if obj.PARAMS.states.RX_DCB && ~(obj.PARAMS.states.RX_DCB_GLO && constIndi == 2) ...
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

if obj.PARAMS.states.RX_DCB_GPS && constIndi == 1 && indGloDcbs(idx)~= 0
    indDcbGlo = obj.INDS_STATE.FLEX_STATES(indGloDcbsi);
    dcbGpsi = x_est_propagated(indDcbGlo);
    Hii(1,indDcbGlo) = 1;
else
    dcbGpsi = 0;
end

if obj.PARAMS.states.RX_DCB_GLO && constIndi == 2 && indGloDcbs(idx)~= 0
    indDcbGlo = obj.INDS_STATE.FLEX_STATES(indGloDcbsi);
    dcbGloi = x_est_propagated(indDcbGlo);
    Hii(1,indDcbGlo) = 1;
else
    dcbGloi = 0;
end

if obj.PARAMS.states.MP_CODE
    indMpCode = obj.INDS_STATE.FLEX_STATES(indMpCodesi);
    mpCodei   = x_est_propagated(indMpCode);
    Hii(1,indMpCode) = 1;
else
    mpCodei = 0;
end

if obj.PARAMS.states.trop
    dtrop = m*state(obj.INDS_STATE.TROP);
    Hii(1,obj.INDS_STATE.TROP) = m;
else
    dtrop = 0;
end

predMeasi = gRange+satBias+rxBias+trop+dtrop+...
    stRangeOffset+relClockCorr+relRangeCorr+ionoCorri+rxDcb+...
    dcbGloi+dcbGpsi+mpCodei;

Hii(1,obj.INDS_STATE.POS)        = A;
Hii(1,obj.INDS_STATE.CLOCK_BIAS(constIndi)) = 1;




sig = 3;









end