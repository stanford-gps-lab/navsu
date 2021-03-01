classdef svOrbitClock < handle

    properties
        % structure containing file directories and ephemeris settings
        settings = struct(...
            'constUse',         [1 0 0 0 0],... % GPS GLO GAL BDS QZSS 
            'miceDir',                   [],... % Directory containing NASA MICE files
            'preciseProdDir',            [],... % Directory containing precise products
            'obsDir',                    [],... % Directory containting IGS observations
            'gpsEphCenter',           'GRG',... % IGS AC code for GPS precise eph
            'gloEphCenter',           'GRG',... % IGS AC code for GLO precise eph
            'galEphCenter',           'GRG',... % IGS AC code for GAL precise eph
            'bdsEphCenter',           'GRG',... % IGS AC code for BDS precise eph
            'gpsClkCenter',           'GRG',... % IGS AC code for GPS precise clk
            'gloClkCenter',           'GRG',... % IGS AC code for GLO precise clk
            'galClkCenter',           'GRG',... % IGS AC code for GAL precise clk
            'bdsClkCenter',           'GRG',... % IGS AC code for BDS precise clk
            'dcbSource',                  2,... % Differential code bias source 2 = DLR- just use that
            'orbitInterpMethod', 'lagrange',...
            'polyfit',   struct('nPolyFit',12,'pfit',8,'cdfit',2),...                 
            'tempDir',                   [],... % directory to put some temporary things...
            'dcbDir',                    [],...; % directory for differential code biases
            'mgxObsDir',                 [],... % MGEX obs files
            'mgxHrObsDir',               []); % MGEX hr obs files
        
        orbMode = 'PRECISE' % whether to use 'PRECISE' or 'BROADCAST' orbit
        clkMode = 'PRECISE' % whether to use 'PRECISE' or 'BROADCAST' clock
        
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
            constUse = res.constUse;
            configFile = res.configFile;
            
            %% Pull info from the .ini file
            if ~isempty(configFile)
                iniData = navsu.thirdparty.ini2struct(configFile);
                basePreciseProdDir = iniData.preciseproddir;
                obsDir             = iniData.obsdir;
            else
                cdi = cd;
                cdi(strfind(cdi,'\')) = '/';
                basePreciseProdDir = [cdi '/data/'];
                obsDir             = [cdi '/data/'];
            end

            obj.settings.preciseProdDir = [basePreciseProdDir 'precise-daily/'];
            obj.settings.mgxObsDir      = [obsDir 'mgex-obs/'];
            obj.settings.mgxHrObsDir    = [obsDir 'mgex-hr-obs/'];
            obj.settings.tempDir        = [basePreciseProdDir 'temp/'];
            obj.settings.dcbDir         = [basePreciseProdDir 'dcb/'];
            obj.settings.navMgxDir      = [basePreciseProdDir 'nav-daily/'];
            
            obj.settings.constUse       = constUse;
            
            %% NASA login/download setup
            defaultCookieFilename = 'nasa_cddis_cookies.txt';
            if ~isempty(res.cookieFile)
                % A specific filename of the cookie was provided
                obj.settings.cookieFile = res.cookieFile;
            elseif ~isempty(res.cookieDir)
                % only a directory was provided-
                obj.settings.cookieFile = fullfile(res.cookieDir,defaultCookieFilename);
            else
                % No location information about where to place teh cookie file was
                % provided
                obj.settings.cookieFile = fullfile(basePreciseProdDir,defaultCookieFilename);
            end
            obj.settings.netrcFile = res.netrcFile;

        end
    end
    
    % function signatures
    methods
        % Load precise orbit and clock
        initOrbitData(obj,varargin)
        
        % Interpolate precise orbit (and clock, but for precise, this is low rate clock)
        [svPos,svVel,iod,svClock,sigma] = propagate(obj,prns,constInds,epochs,varargin)
        
        % Load precise clock (higher rate than clock with orbits)
        initClockData(obj,year,doy,varargin)
        
        % Interpolate high rate precise clock
        cbias = clock(obj,prns,constInds,epochs);
        
        % Load IGS TEC map
        initIonoData(obj,year,doy,varargin);
        
        % Iono delay from TEC map
        [tecs,delays,tecSlant] = ionoDelay(obj,epoch,llh,varargin)    
        
        % Load atx file
        initAtxData(obj,filenameAtx);
        
        % Load DCB data
        initDcb(obj,year,doy)
        
        % Use the clock data from the .sp3 to populate the precise clock
        % field
        initPClockFromPEph(obj)
        
        % Load broadcast data
        initBroadcastData(obj,years,doys,varargin)
        
        %% The rest of the functions mostly help with the above functions. 
        % Helps with precise clock interpolation
        cbias = clockInterp(obj,prns,constInds,epochs,Clck);
        % Helps with precise clock interpolation
        cbias = clockBiasFromProd(obj,prns,constInds,epochs)
        
        % Helps with precise orbit interpolation
        [svPos,svVel,iod,sigOrbit] = svPosFromProd(obj,prns, epochs,settings,...
            pPosInds,pPosPoly,constInds,FLAG_APC_OFFSET,atxData,sunPos,dttx);
        % Helps with precise orbit interpolation
        [posP, velP,pPosInds,pPosPoly,ROut] = PPosInterp(obj,PRNs, epochs,Ppos,Pprns,...
            Pepochs,settings,Pvel,pPosInds,pPosPoly,constInds,PconstInds,FLAG_APC_OFFSET,...
            atxData,sunPos,dttx)
        
        [svPos,svVel,iod,svClock] = predictOrbit(obj,prns,constInds,epochs,latency);
        [yy, yy_dot, chi2, p] = polyinterp(obj,x, y, m_order, xx, flag, var,polyIn)
        cbias = predictClock(obj,prns,constInds,epochs,latency);
        
    end


end















