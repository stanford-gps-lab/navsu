classdef pppFilter < matlab.mixin.Copyable
    
    
    properties
        
        state   % the state of the filter -> should be of type pppanal.ppp.State
        cov     % covariance of the state
        
        pos                % ECEF position
        vel                % ECEF velocity
        R_b_e              % DCM from body to ECEF
        imuBiasStates      % imu bias states
        clockBias          % receiver clock bias(es)
        clockDrift         % receiever clock drift(s)
        carrierAmbiguities % carrier phase ambiguity estimates
        
        INDS_STATE % state indexing
        
        lastAccMeas % most recent acceleration measurements
        lastGyroMeas % most recent gyro measurement
        epochLastInertialUpdate % time of the last inertial update
        epochLastGnssUpdate      % time of the last GNSS update
        
        posPrevTc % position at previous tight coupling update
        
        % contains latest geometry free combinations
        cycleSlipInfo = struct('gFree',zeros(0,1),'epochLastGFree',zeros(0,1),'measInfoGFree',zeros(0,4));
        
        % phase windup
        phWind = struct('phaseOffset',zeros(0,1),'PrnConstInd',zeros(0,2));
        
        StateMap % mapping matrix- indicates what the state is in each position of the state and covariance matrix
        
        % all satellites used in the solution- useful for solution
        % separation :)
        allSatsSeen
        
        initialized = false % whether or not the fitler has been initialized
        
        PARAMS % parameters associated with the running of this filter!
        
        resids % extra info for output about measurement residuals
        
        measRemoved % extra info for measurements that were removed :)
        
    end
    
    
    methods
        function obj = pppFilter()
            % TO DO: ALLOW FOR THE SETTING OF THE PARAMETERS WITH A
            % CONFIGURATION FILE OR SOMETHING
            
            obj.PARAMS = obj.initParams;
            
        end
    end
    
    
    % function signatures
    methods
        complete = initialize(obj,corrData,obs,varargin)
        
        % the time AND measurement update :O
        [measMatRemoved,measMatRemovedLow] = update(obj,epoch,obs,corrData)
        
        manageStates(obj,epoch,gnssMeas,PARAMS,outStruc);
        measRemoved = checkCycleSlips(obj,epoch,gnssMeas,PARAMS);
        removeFlexState(obj,measInfoRemove);
        
        % This function basically just stores the default configuration of
        % the filter
        PARAMS = initParams(obj);
        
        outData = saveState(obj,outData,epoch,obs);
        
    end
    
    methods (Access = private)
        timeUpdate(obj,epoch)
        
        measUpdate(obj,epoch,obs,corrData,measRemovedSlip)
        
        [predMeas,H,sig] = doppModel(obj,nState,dVel,A,rxDrift,constInd)
        
        [predMeas,H,sig] = carrierModel(obj,nState,sigi,freqi,tecSlant,state,m,indIonosi, ...
            indMpCarrsi,indAmbStatesi,phWind,gRange,satBias,rxBias,trop,stRangeOffset,...
            relClockCorr,relRangeCorr,A,constIndi)
        
        [predMeas,H,sig] = codeModel(obj,nState,sigi,freqi,tecSlant,state,...
            constIndi,indGloDcbsi,indMpCodesi,m,gRange,satBias,rxBias,trop,stRangeOffset,...
            relClockCorr,relRangeCorr,A)
        
        [pred_meas,H,R,el,az,prnConstInds,measMatRemovedLow,measMat] = handleGnssMeas(obj,epoch,obs,corrData)
       
        [pred_measi,Hi,ri,measMati] = handleVehicleConstraintPseudomeas(obj)
        
        [predMeasi,Hi,Ri,measMati] = handlePositionMeas(obj,posMeas)

    end
    
    methods(Static)
        plotOutput(outputs,varargin)
    end
    
    
end















