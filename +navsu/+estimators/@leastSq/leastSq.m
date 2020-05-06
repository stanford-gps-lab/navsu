classdef leastSq < matlab.mixin.Copyable
    
    
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
        
        % state indices
        INDS_STATE
        
        % all satellites used in the solution- useful for solution
        % separation :)
        allSatsSeen
                
        initialized = true % it's always initialized :)
        
        PARAMS % parameters associated with the running of this filter!
        
        resids % extra info for output about measurement residuals
        
        measRemoved % extra info for measurements that were removed :)
        
    end
    
    
    methods
        function obj = leastSq()
            % TO DO: ALLOW FOR THE SETTING OF THE PARAMETERS WITH A
            % CONFIGURATION FILE OR SOMETHING
            
            obj.PARAMS = obj.initParams;
        end
    end
    
    
    % function signatures
    methods
        complete = initialize(obj,corrData,varargin)
        
        % the time AND measurement update :O
        update(obj,epoch,obs,corrData)
        
        % This function basically just stores the default configuration of
        % the filter
        PARAMS = initParams(obj);
        
        outData = saveState(obj,outData,epoch,obs);
        
        obs = checkMeas(obj,obs0)
    end
    
    methods(Static)
        plotOutput(outputs,varargin)
    end
    
    
end















