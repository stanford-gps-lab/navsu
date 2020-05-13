classdef leastSq < navsu.estimators.AbstractNavFilter
    
    
    properties        
        % all satellites used in the solution- useful for solution
        % separation :)
        allSatsSeen
                
%         initialized = true % it's always initialized :)
        
%         PARAMS % parameters associated with the running of this filter!
        
        resids % extra info for output about measurement residuals
        
        measRemoved % extra info for measurements that were removed :)
        
    end
    
    
    methods
        function obj = leastSq()
            % TO DO: ALLOW FOR THE SETTING OF THE PARAMETERS WITH A
            % CONFIGURATION FILE OR SOMETHING
            
            obj.PARAMS = obj.initParams;
            
            obj.initialized = true;
        end
    end
    
    
    % function signatures
    methods
        complete = initialize(obj,obs,corrData,varargin)
        
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















