classdef leastSq < navsu.estimators.pppFilter
    
       
    
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
        measId = initialize(obj,obs,corrData,varargin)
        
        % the time AND measurement update :O
        [measId,extraInputs] = update(obj,epoch,obs,corrData,varargin)
        
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















