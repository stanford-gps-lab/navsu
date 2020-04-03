classdef svOrbitClock < handle


    properties

        orbMode = 'PRECISE' % whether to use precise or broadcast orbit
        clkMode = 'PRECISE' % whether to use precise or broadcast clock
        
        PClock % precise clock data
        PEph   % precise orbit data
        
        BEph   % broadcast orbit and clock data
        
        iono   % ionospheric data- could be TEC map
        
        atx    % Antenna phase center file from IGS- parased
        
        dcb    % differential code bias information
        
        settings = initSettings% its the full initSettings structure
    end


    methods
        function obj = svOrbitClock(varargin)
            
            p = inputParser;
            
            p.addParameter('settings',initSettings);
            
            % parse the results
            parse(p, varargin{:});
            res = p.Results;
            settings = res.settings;
            
            obj.settings = settings;
        end
    end
    
    % function signatures
    methods
        % Load precise orbit and clock
        initPEph(obj,varargin)
        % Interpolate precise orbit
        [svPos,svVel,iod,sigEph] = propagate(obj,prns,constInds,epochs,varargin)
        % Load precise clock
        initPClock(obj,year,doy,varargin)
        % Interpolate precise clock
        cbias = clock(obj,prns,constInds,epochs);
        % Load IGS TEC map
        initIonoData(obj,year,doy,varargin);
        % Iono delay from TEC map
        [tecs,delays,tecSlant] = ionoDelay(obj,epoch,llh,varargin)    
        % Load atx file
        initAtxData(obj,filenameAtx);
        % Load differential code bias data
        initDcb(obj,year,doy);
    end


end















