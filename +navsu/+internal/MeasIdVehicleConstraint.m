classdef MeasIdVehicleConstraint < navsu.internal.MeasID
    properties (SetObservable = true )       
        id = 255;
        
        subtype navsu.internal.MeasEnum;
    end
    
    methods
        function obj = MeasIdVehicleConstraint(id,subtype)
            % Instantiate the array of the super class
            obj = obj@navsu.internal.MeasID(repelem(navsu.internal.MeasEnum.VehicleConstraint,...
                size(id,1),size(id,2)));
            for idx = 1:size(id,1)
                for jdx = 1:size(id,2)
                    
                    obj(idx,jdx).TypeID = navsu.internal.MeasEnum.VehicleConstraint;
                    
                    obj(idx,jdx).id = id(idx,jdx);
                    
                    obj(idx,jdx).subtype = subtype(idx,jdx);
                    
                    obj(idx,jdx).idVec = permute([double(navsu.internal.MeasEnum.VehicleConstraint) ...
                        id(idx,jdx) double(subtype(idx,jdx)) 0 0 0],[1 3 2]);       
                    
                end
            end
        end
    end
    
    
end
