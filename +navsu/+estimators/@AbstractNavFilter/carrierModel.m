function [predMeas,H,sig] = carrierModel(obj,nState,sigi,freqi,tecSlant,state,m,indIonosi, ...
    indMpCarrsi,indAmbStatesi,phWind,gRange,satBias,rxBias,trop,stRangeOffset,...
    relClockCorr,relRangeCorr,A,constIndi)
predMeas = [];
H = [];
sig= [];

H = zeros(1,nState);

%% carrier
if sigi < 100 && strcmp(obj.PARAMS.states.ionoMode,'TEC')
    % If it's single frequency, need to include iono correction
    ionoCorri = +tecSlant*40.3*10^15./freqi.^2;
    
elseif sigi < 100 && strcmp(obj.PARAMS.states.ionoMode,'L1DELAYSTATE')
    ionoCorrModel = -tecSlant*40.3*10^15./freqi.^2;
    
    indIono = obj.INDS_STATE.FLEX_STATES(indIonosi);
    delayL1i = state(indIono);
    
    hi = (1575.42e6).^2./freqi.^2;
    H(1,indIono) = hi;
    
    ionoCorri = hi*delayL1i+ionoCorrModel;
else
    ionoCorri = 0;
end

if obj.PARAMS.states.MP_CARR
    indMpCarr = obj.INDS_STATE.FLEX_STATES(indMpCarrsi);
    mpCarri   = state(indMpCarr);
    H(1,indMpCarr) = 1;
else
    mpCarri = 0;
end

if obj.PARAMS.states.trop
    dtrop = m*state(obj.INDS_STATE.TROP);
    H(1,obj.INDS_STATE.TROP) = m;
else
    dtrop = 0;
end


indAmbState = obj.INDS_STATE.FLEX_STATES(indAmbStatesi);

ambEst = state(indAmbState);

% Carrier phase windup
phWindi = phWind*navsu.constants.c/freqi;

predMeas = gRange+satBias+rxBias+...
    trop+dtrop+stRangeOffset+relClockCorr+...
    relRangeCorr+ionoCorri+ambEst+phWindi+mpCarri;

H(1,obj.INDS_STATE.POS)        = A;
H(1,obj.INDS_STATE.CLOCK_BIAS(constIndi)) = 1;
H(1,indAmbState)              = 1;

% Also include lever arm sensitivity lol
% just constructing things

% z axis (
rArmBody = obj.PARAMS.IMU_ARM;
if norm(rArmBody) > 0
    
    rArmBodyNorm = rArmBody./norm(rArmBody);
    
    % Convert the LOS to the body frame
    rhoBody = obj.R_b_e'*A';
    
    Hangle(1) = sum(norm(rArmBody([2 3])).*cross([1 0 0]',rArmBodyNorm).*rhoBody);
    Hangle(2) = sum(norm(rArmBody([1 3])).*cross([0 1 0]',rArmBodyNorm).*rhoBody);
    Hangle(3) = sum(norm(rArmBody([1 2])).*cross([0 0 1]',rArmBodyNorm).*rhoBody);
    H(1,obj.INDS_STATE.ATTITUDE) = Hangle;
    
    % ri = weighti.*0.03^2;
    
end
sig = 0.003;













end