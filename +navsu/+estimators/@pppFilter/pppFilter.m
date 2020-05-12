classdef pppFilter < navsu.estimators.AbstractNavFilter
    
    
    properties
        
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
    
    methods %(Access = private)
        timeUpdate(obj,epoch)
        
        measUpdate(obj,epoch,obs,corrData,measRemovedSlip)
        
     
%         [posApc,velApc] = posVelApc(obj);  % Position and velocity of the GNSS antenna phase center. 
    end
    
    methods(Static)
        plotOutput(outputs,varargin)
    end
    
    
end















