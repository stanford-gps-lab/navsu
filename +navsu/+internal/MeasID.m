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
    
    methods (Sealed)
        function matchVec = matches(obj,measIdList)
            
            if ~iscolumn(measIdList)
                error('Input must be column vector')
            end
            
            if length(obj) > 1
                error('Can only check one at at time right now')
            end
            
            matchVec = zeros(size(measIdList));
            
            types = cat(1,measIdList.TypeID);
            indsTypeMatch = find(types == obj.TypeID);
            
            % Only check measurements of the same type from now on
            if isempty(indsTypeMatch)
                return;
            end
               
            measSameType = measIdList(indsTypeMatch);
            
            matchVecSame = ones(length(indsTypeMatch),1);
            
            switch obj.TypeID
                case navsu.internal.MeasEnum.GNSS
                     prnList = cat(1,measSameType.prn);
                     constList = cat(1,measSameType.const);
                     freqList  = cat(1,measSameType.freq);
                     subtypeList = cat(1,measSameType.subtype);
                     
                     % prn check
                     matchVecSame = matchVecSame & ((obj.prn == 0) | (prnList == 0) | ...
                         obj.prn == prnList);
                     
                     % constellation check
                     matchVecSame = matchVecSame & ((obj.const == 0) | (constList == 0) | ...
                         obj.const == constList);
                     
                      % freq check
                     matchVecSame = matchVecSame & ((obj.freq == 0) | (freqList == 0) | ...
                         obj.freq == freqList);
                    
                     % subtype check
                     matchVecSame = matchVecSame & ((obj.subtype == 0) | (subtypeList == 0) | ...
                         obj.subtype == subtypeList);
                    
                case navsu.internal.MeasEnum.Position
                    xyzList = cat(1,measSameType.xyz);
                    idList = cat(1,measSameType.id);
                    
                    % xyz check
                    matchVecSame = matchVecSame & ((obj.xyz == 0) | (xyzList == 0) | ...
                        obj.xyz == xyzList);
                    
                    % id check
                    matchVecSame = matchVecSame & ((obj.id == 0) | (idList == 0) | ...
                        obj.id == idList);
                                        
                case navsu.internal.MeasEnum.Velocity
                    xyzList = cat(1,measSameType.xyz);
                    idList = cat(1,measSameType.id);
                    
                    % xyz check
                    matchVecSame = matchVecSame & ((obj.xyz == 0) | (xyzList == 0) | ...
                        obj.xyz == xyzList);
                    
                    % id check
                    matchVecSame = matchVecSame & ((obj.id == 0) | (idList == 0) | ...
                        obj.id == idList);
                    
                case navsu.internal.MeasEnum.IMU
                    error('implement me')
                    
                case navsu.internal.MeasEnum.VehicleConstraint
                    error('implement me')
                    
            end
            
            matchVec(indsTypeMatch) = matchVecSame;
            
        end
    end
    
    
    
end
    
