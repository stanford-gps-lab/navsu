classdef svOrbitClock < handle

    properties
        % structure containing file directories and ephemeris settings
        settings % see navsu.internal.initSettings
        
        orbMode(1,:) char {mustBeMember(orbMode, {'PRECISE', 'BROADCAST', 'PREDICT'})} = 'PRECISE' % 'PREDICT" is deprecated
        clkMode(1,:) char {mustBeMember(clkMode, {'PRECISE', 'BROADCAST', 'PREDICT'})} = 'PRECISE' % 'PREDICT" is deprecated
        
        PClock % precise clock data structure
        PEph   % precise orbit data structure
        
        BEph   % broadcast orbit and clock data 
        
        iono   % ionospheric data- could be TEC map
        
        atx    % Antenna phase center file from IGS- parsed
        
        dcb    % GNSS satellite differential code bias data
    end

    methods
        function obj = svOrbitClock(varargin)
            % svOrbitClock
            % DESCRIPTION:
            %   Class to handle products for GNSS solutions.  This includes
            %   orbit and clock (precise and/or broadcast), ionospheric
            %   maps, antenna phase center data, and differential code bias
            %   data.  This includes functionality to download the data as
            %   well as use them to apply corrections. 
            %   Also contains bunch of information about local product file
            %   structure. 
            %   
            % OPTIONAL INPUT:
            %   constUse   - Constellation usage mask 1x5 boolean vector,
            %                with GRECS for each spot.  For example, 
            %                [1 1 0 0 0] means GPS and GLONASS are enabled
            %   configFile - configuration file.  Should match the file
            %                default.ini file, but please name it something
            %                else. This tells the code where to put your
            %                products on your local machine when they are
            %                downloaded.
            %   netrcFile  - REQUIRED FOR NASA CDDIS DOWNLOADS. 
            %                Username and password file for NASA CDDIS 
            %                data/product access. Information about general 
            %                access can be found at [1], and specifics of
            %                what this file should look like can be found at 
            %                [2] once a username and password have been created.
            %   cookieDir  - Local folder where cookies from NASA CDDIS web 
            %                access can/should be stored.
            %   cookieFile - A storage file can also be specified for cookie storage. 
            % 
            % References: 
            % [1] https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
            % [2] https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html
            %
            % OUTPUT:
            %   obj        - the object!
            %
            % See also:  the rest of the methods in here!
                      
            %% Parse inputs
            p = inputParser;
            p.addParameter('constUse',[1 0 0 0 0]);
            p.addParameter('configFile',[]);
            p.addParameter('netrcFile',[]);
            p.addParameter('cookieFile',[]);
            p.addParameter('cookieDir',[]);
            parse(p, varargin{:});
            res = p.Results;
            
            %% initialize settings
            obj.settings = navsu.internal.initSettings('configFile', res.configFile, ...
                                                       'netrcFile', res.netrcFile, ...
                                                       'cookieFile', res.cookieFile, ...
                                                       'cookieDir', res.cookieDir);
            obj.settings.constUse = res.constUse;

            % a few fields now need to be overwritten to be backwards
            % compatible:

            % remove the '/mgex/' part from the path
            obj.settings.navMgxDir = fullfile(obj.settings.baseDir, 'nav-daily/');
            
            % reset the ephemeris and clk center
            centerConsts = {'gps'; 'glo'; 'gal'; 'bds'};
            for cci = 1:length(centerConsts)
                obj.settings.([centerConsts{cci}, 'EphCenter']) = 'GRG';
                obj.settings.([centerConsts{cci}, 'ClkCenter']) = 'GRG';
            end

        end
    end
    
    % function signatures
    methods
        %% Initialization methods
        % Load precise orbit and clock
        initOrbitData(obj, varargin)
        
        % Load precise clock (higher rate than clock with orbits)
        initClockData(obj, year, doy, varargin)
        
        % Load IGS TEC map
        initIonoData(obj, year, doy, varargin);
        
        % Iono delay from TEC map
        [tecs, delays, tecSlant] = ionoDelay(obj, epoch, llh, varargin)    
        
        % Load atx file
        initAtxData(obj, filenameAtx);
        
        % Load DCB data
        initDcb(obj, year, doy)
        
        % Use the clock data from the .sp3 to populate the precise clock
        % field
        initPClockFromPEph(obj)
        
        % Load broadcast data
        initBroadcastData(obj, years, doys, varargin)
        
        %% Main propagation methods
        % Interpolate precise orbit (and clock, but for precise, this is 
        % low rate clock)
        [svPos, svVel, iod, svClock, sigma] = propagate(obj, ...
            prns, constInds, epochs,varargin)
        
        % Interpolate high rate precise clock
        cbias = clock(obj, prns, constInds, epochs);
        
        %% The rest of the functions mostly help with the above functions. 
        % Helps with precise clock interpolation
        cbias = clockInterp(obj, prns, constInds, epochs);
        
        % Helps with precise orbit interpolation
        [posP, velP, pPosInds, pPosPoly, ROut] = PPosInterp(obj, prns, ...
            constInds, epochs, pPosInds, pPosPoly, FLAG_APC_OFFSET, ...
            sunPos, dttx);
        
        [svPos, svVel, iod, svClock] = predictOrbit(obj, prns, ...
            constInds, epochs, latency);
        
        cbias = predictClock(obj, prns, constInds, epochs, latency);
                
    end

end















