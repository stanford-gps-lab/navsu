classdef MeasIdPos < navsu.internal.MeasID
    
    
    properties (SetObservable = true)
        
        id = 255;
        
        xyz = 255;
        
    end
    
    
    methods
        function obj = MeasIdPos(id,xyz)
            
            % Instantiate the array of the super class
            obj = obj@navsu.internal.MeasID(repelem(navsu.internal.MeasEnum.Position,...
                size(id,1),size(id,2)));
            for idx = 1:size(id,1)
                for jdx = 1:size(id,2)
                    
                    obj(idx,jdx).TypeID = navsu.internal.MeasEnum.Position;
                    
                    obj(idx,jdx).id = id(idx,jdx);
                    
                    obj(idx,jdx).xyz = xyz(idx,jdx);
                end
            end
        end
    end
    
    
end









