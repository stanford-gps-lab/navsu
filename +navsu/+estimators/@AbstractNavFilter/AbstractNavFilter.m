classdef (Abstract) AbstractNavFilter < matlab.mixin.Copyable 
    
    properties
        initialized = 0 % whether or not the filter has been initialized.
        % 0 = not started. 1 = initialized but that's it. 2 = initialized
        % with at least one time and measurement update
        
    end
    
    methods (Abstract)
        % These MUST be implemented! 
        obs = checkMeas(obj,obs);
        
        % Initialization
        measId = initialize(obj,obs,corrData,varargin)
        
        % the time AND measurement update :O
        [measId,measMatRemoved,measMatRemovedLow] = update(obj,epoch,obs,corrData,varargin)
        
        % Save what you would like 
        outData = saveState(obj,outData,epoch,obs);
    end    
    
    
end
    