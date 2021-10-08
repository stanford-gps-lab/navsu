classdef pppFilter < navsu.estimators.AbstractNavFilter
    
    
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
        
        
        PARAMS % parameters associated with the running of this filter!
        
         % phase windup
        phWind = struct('phaseOffset',zeros(0,1),'PrnConstInd',zeros(0,2));
        
        
        lastAccMeas % most recent acceleration measurements
        lastGyroMeas % most recent gyro measurement
        epochLastInertialUpdate % time of the last inertial update
        epochLastGnssUpdate      % time of the last GNSS update
        
        posPrevTc % position at previous tight coupling update
        
        % contains latest geometry free combinations
        cycleSlipInfo = struct('gFree',zeros(0,1),'epochLastGFree',zeros(0,1),'measInfoGFree',zeros(0,4));
                
        StateMap % mapping matrix- indicates what the state is in each position of the state and covariance matrix
        
        % all satellites used in the solution- useful for solution
        % separation :)
        allSatsSeen
           
        resids % extra info for output about measurement residuals
        
        measRemoved % extra info for measurements that were removed :)
    end
    
    properties (Constant)
        % Reasons for removing a measurement
        RemovedLow = 1;
        RemovedResid = 2;
        RemovedSlip  = 3;
        RemovedEl    = 4;
        RemovedCn0   = 5;
        
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
        
        obs = checkMeas(obj,obs);
        
        measId = initialize(obj,corrData,obs,varargin)
        
        % the time AND measurement update :O
        [measId,extraInputs] = update(obj,epoch,obs,corrData,varargin)
        
        measRemovedSlip = manageStates(obj,epoch,gnssMeas,PARAMS,outStruc)
        measRemovedSlip = manageStatesMulti(obj,epoch,obs);
        measRemoved = checkCycleSlips(obj,epoch,gnssMeas,PARAMS);
        removeFlexState(obj,measInfoRemove);
        
        % This function basically just stores the default configuration of
        % the filter
        PARAMS = initParams(obj);
        
        outData = saveState(obj,outData,epoch,obs);
               
    end
    
    
    methods
        % These classes have already been implemented and are mostly just
        % useful for the various GNSS estimators, including least squares
        % and inertial PPP filters
        
        % Build predicted measurements, sensitivity matrix, and pull
        % measurement IDs for GNSS measurements
        [predMeas,H,R,el,az,prnConstInds,idList,measList,measIdRemovedLow,...
            extraInputs,dop] =  handleGnssMeas(obj,epoch,obs,corrData,varargin)
        
        % doppModel- called within handleGnssMeas and produces predicted
        % value and sensitivity matrix for doppler measurements
        [predMeas,H,sig] = doppModel(obj,nState,dVel,A,rxDrift,constInd)
        
        % carrierModel- called within handleGnssMeas and produces predicted
        % value and sensitivity matrix for carrier phase measurements
        [predMeas,H,sig] = carrierModel(obj,nState,sigi,freqi,tecSlant,state,m,indIonosi, ...
            indMpCarrsi,indAmbStatesi,phWind,gRange,satBias,rxBias,trop,stRangeOffset,...
            relClockCorr,relRangeCorr,A,constIndi,indEphErri)
        
        % codeModel- called within handleGnssMeas and produces predicted
        % value and sensitivity matrix for code phase measurements
        [predMeas,H,sig] = codeModel(obj,SimpleModel,nState,sigi,freqi,tecSlant,state,...
            constIndi,indGloDcbsi,indMpCodesi,m,gRange,satBias,rxBias,trop,stRangeOffset,...
            relClockCorr,relRangeCorr,A,indIonosi,indEphErri)
        
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
        [complete, measId,extraInputs] = leastSquaresSol(obj,epoch,obs,corrData,varargin)

        [posApc,velApc] = posVelApc(obj);  % Position and velocity of the GNSS antenna phase center.

        outState = saveOutStatePpp(obj,outData,epoch,obs)
        
        plotOutput(bbj,outputs,varargin)
        
        plotOutputPpp(obj,outputs,varargin)
    end
    
    methods %(Access = private)
        timeUpdate(obj,epoch)
        
        [measId,extraInputs] = measUpdate(obj,epoch,obs,corrData,measRemovedSlip,varargin)
        
    end
    
end















