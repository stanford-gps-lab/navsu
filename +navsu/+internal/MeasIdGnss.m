classdef MeasIdGnss < navsu.internal.MeasID
    
    % GNSS Measurement ID
    properties (SetObservable = true, SetAccess = protected)
        
        prn    = 0;
        
        const  = 0;
        
        freq   = 0;
        
        subtype navsu.internal.MeasEnum
                
    end
    
    
    methods
        function obj = MeasIdGnss(prn,const,freq,subtype)

            % Instantiate the array of the super class 
            obj = obj@navsu.internal.MeasID(repelem(navsu.internal.MeasEnum.GNSS,...
                size(prn,1),size(prn,2)));
            for idx = 1:size(prn,1)
                for jdx = 1:size(prn,2)
                    
                    obj(idx,jdx).TypeID = navsu.internal.MeasEnum.GNSS;  % 8 bits
                    
                    obj(idx,jdx).prn = prn(idx,jdx);   % 
                    
                    obj(idx,jdx).const   = const(idx,jdx);
                    
                    obj(idx,jdx).freq    = freq(idx,jdx);
                    
                    obj(idx,jdx).subtype = subtype(idx,jdx);
                    
                    obj(idx,jdx).idVec = permute([double(navsu.internal.MeasEnum.GNSS) ...
                        prn(idx,jdx) const(idx,jdx) freq(idx,jdx) double(subtype(idx,jdx)) 0],[1 3 2]);
                end
            end
            
            
            
        end
        
        
        
    end
    
    
    
    
    
    
    
    
    
end