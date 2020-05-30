classdef MeasIdWheels < navsu.internal.MeasID
    properties (SetObservable = true )       
        id = 255;
        
        subtype navsu.internal.MeasEnum;
    end
    
    methods
        function obj = MeasIdWheels(id,subtype)
            % Instantiate the array of the super class
            obj = obj@navsu.internal.MeasID(repelem(navsu.internal.MeasEnum.Wheels,...
                size(id,1),size(id,2)));
            for idx = 1:size(id,1)
                for jdx = 1:size(id,2)
                    obj(idx,jdx).id = id(idx,jdx);
                    
                    obj(idx,jdx).subtype = subtype(idx,jdx);
                end
            end
        end
    end
    
end









