classdef MeasIdPos < navsu.internal.MeasID
    
    
    properties (SetObservable = true)
               
       id
       
       xyz
        
    end
    
    
    methods 
        function obj = MeasIdPos(id,xyz)
            obj.TypeID = navsu.internal.MeasEnum.Position;
            
            obj.id = id;
            
            obj.xyz = xyz;
            
        end
        
        
    end
    
    
    
    
    
    
    
    
    
end