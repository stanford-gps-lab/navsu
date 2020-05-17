classdef inertialPppFilter < navsu.estimators.pppFilter
    
    properties
        
        
    end
    
    methods
         % the time AND measurement update :O
        [measMatRemoved,measMatRemovedLow] = update(obj,epoch,obs,corrData,varargin)
        
        % IMU mechanization
        mechanization(obj,epoch,obs);
        
        % time update is different from standard ppp
        timeUpdate(obj,epoch,obs)
        
        outData = saveState(obj,outData,epoch,obs)
        
    end
    
    
    methods(Static)
        plotOutputInertial(outputs,varargin)
    end
    
end















