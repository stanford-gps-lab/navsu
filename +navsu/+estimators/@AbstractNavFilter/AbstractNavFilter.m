classdef (Abstract) AbstractNavFilter < matlab.mixin.Copyable 
    
    properties
        state   % the state of the filter 
        cov     % covariance of the state
        
        pos                % ECEF position
        vel                % ECEF velocity
        R_b_e              % DCM from body to ECEF
        imuBiasStates      % imu bias states
        clockBias          % receiver clock bias(es)
        clockDrift         % receiever clock drift(s)
        carrierAmbiguities % carrier phase ambiguity estimates
        
        INDS_STATE % state indexing
        
        initialized = false % whether or not the fitler has been initialized
        
        PARAMS % parameters associated with the running of this filter!
        
         % phase windup
        phWind = struct('phaseOffset',zeros(0,1),'PrnConstInd',zeros(0,2));
        
    end
    
    methods (Abstract)
        % These MUST be implemented! 
        obs = checkMeas(obj,obs);
        
        % Initialization
        complete = initialize(obj,corrData,obs,varargin)
        
        % the time AND measurement update :O
        [measMatRemoved,measMatRemovedLow] = update(obj,epoch,obs,corrData)
        
        % Save what you would like 
        outData = saveState(obj,outData,epoch,obs);
    end    
    
    methods
        % These classes have already been implemented and are mostly just
        % useful for the various GNSS estimators
        
        
        % Build predicted measurements, sensitivity matrix, and pull
        % measurement IDs for GNSS measurements
        [predMeas,H,R,el,az,prnConstInds,idList,measList,measIdRemovedLow] =  ...
            handleGnssMeas(obj,epoch,obs,corrData,varargin)
        
        % doppModel- called within handleGnssMeas and produces predicted
        % value and sensitivity matrix for doppler measurements
        [predMeas,H,sig] = doppModel(obj,nState,dVel,A,rxDrift,constInd)
        
        % carrierModel- called within handleGnssMeas and produces predicted
        % value and sensitivity matrix for carrier phase measurements
        [predMeas,H,sig] = carrierModel(obj,nState,sigi,freqi,tecSlant,state,m,indIonosi, ...
            indMpCarrsi,indAmbStatesi,phWind,gRange,satBias,rxBias,trop,stRangeOffset,...
            relClockCorr,relRangeCorr,A,constIndi)
        
        % codeModel- called within handleGnssMeas and produces predicted
        % value and sensitivity matrix for code phase measurements
        [predMeas,H,sig] = codeModel(obj,SimpleModel,nState,sigi,freqi,tecSlant,state,...
            constIndi,indGloDcbsi,indMpCodesi,m,gRange,satBias,rxBias,trop,stRangeOffset,...
            relClockCorr,relRangeCorr,A)
        
        % handleVehicleConstraintPseudomeas- produce predicted value and
        % sensitivity matrix for vehicle slip constraints
        [pred_measi,Hi,ri,measIdi,measi] = handleVehicleConstraintPseudomeas(obj)
        
        % handlePositionMeas- produce predicted value and
        % sensitivity matrix for direct position measurements
        [predMeasi,Hi,Ri,measIdi,measi] = handlePositionMeas(obj,posMeas)
        
        % handleVelocityMeas- produce predicted value and
        % sensitivity matrix for direct velocity measurements
        [predMeasi,Hi,Ri,measIdi,measi] = handleVelocityMeas(obj,velMeas)
        
        % leastSquaresSol - produce a least squares solution based on the 
        % available observations.  This currently requires GNSS
        % observations
        [complete, measId] = leastSquaresSol(obj,epoch,obs,corrData)

        [posApc,velApc] = posVelApc(obj);  % Position and velocity of the GNSS antenna phase center.

    end
    
end
    