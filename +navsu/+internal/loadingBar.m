classdef loadingBar < handle
    
    properties
        
        pctDone = 0 % integer percent completed
        runTimeStart % 'tic' at the initialization
        
        figHandle % the figure handle
        
        nTimes % number of total times
        
        pctUpdate = 1 % update every this percentage
        
        DisplayFlag = true;
    end
    
    
    methods
        function obj = loadingBar(nTimes)
            obj.nTimes = nTimes;
            
            obj.runTimeStart = tic;
            
            obj.figHandle = waitbar(0,'0 Percent Complete');
            
        end
        
        function update(obj,tIndex)
            if mod(floor(tIndex/obj.nTimes*100),obj.pctUpdate) == 0 && ...
                    floor(tIndex/obj.nTimes*100) > obj.pctDone
                tElapsed = toc(obj.runTimeStart);
                tRemaining = tElapsed*(obj.nTimes-tIndex)./tIndex;
                obj.pctDone =  floor(tIndex/obj.nTimes*100);
                waitbar(obj.pctDone/100,obj.figHandle,[num2str(obj.pctDone) '% Complete, '...
                    num2str(tRemaining/60,'%5.2f') ' Minutes Remaining']);
            end
            
        end
        
        function close(obj);
           close(obj.figHandle) 
        end
    end
    
    
    
    
end







