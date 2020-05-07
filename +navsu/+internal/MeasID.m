classdef  MeasID < matlab.mixin.Heterogeneous
    
    % Heterogenous measurement ID super class
    properties
        
        % Type of measurement
        TypeID navsu.internal.MeasEnum
        
        
    end
    
    methods
        
        function obj = MeasID(type)
            % Initial constructor
            if nargin > 0
                obj = repelem(obj,size(type,1),size(type,2));
                for idx = 1:size(type,1)
                    for jdx = 1:size(type,2)
                        obj(idx,jdx).TypeID = type(idx,jdx);
                    end
                end
            else
                obj.TypeID = [];
            end
            
        end
     
    end
    
end
    
