classdef MeasIdGnss < navsu.internal.MeasID
    
    % GNSS Measurement ID
    properties (SetObservable = true)
        
        prn    = 255;
        
        const  = 255;
        
        freq   = 255;
        
        subtype = 255;
                
    end
    
    
    methods
        function obj = MeasIdGnss(prn,const,freq,subtype)

            % Instantiate the array of the super class 
            obj = obj@navsu.internal.MeasID(repelem(navsu.internal.MeasEnum.GNSS,...
                size(prn,1),size(prn,2)));
            for idx = 1:size(prn,1)
                for jdx = 1:size(prn,2)
                    
                    obj(idx,jdx).TypeID = navsu.internal.MeasEnum.GNSS;
                    
                    obj(idx,jdx).prn = prn(idx,jdx);
                    
                    obj(idx,jdx).const   = const(idx,jdx);
                    
                    obj(idx,jdx).freq    = freq(idx,jdx);
                    
                    obj(idx,jdx).subtype = subtype(idx,jdx);

                end
            end
            
        end
        
        
        
    end
    
    
    
    
    
    
    
    
    
end