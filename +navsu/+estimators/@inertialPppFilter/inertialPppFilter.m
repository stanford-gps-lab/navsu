classdef inertialPppFilter < navsu.estimators.pppFilter
    
    properties
        
        
    end
    
    methods
         % the time AND measurement update :O
        [measMatRemoved,measMatRemovedLow] = update(obj,epoch,obs,corrData)
        
        % IMU mechanization
        mechanization(obj,epoch,obs);
        
        % time update is different from standard ppp
        timeUpdate(obj,epoch,obs)
        
        outData = saveState(obj,outData,epoch,obs)
        
    end
    
end















