function [dayChangei,YearChangei] = download(inChoice,YearList,dayList,settings,varargin)
% download
% Downloads various IGS data and products from FTP sites. 
%
% INPUTS:
%  inChoice            - selection of what to download- options below
%  YearList            - N-length vector indicating desired years of
%                        desired download
%  dayList             - N-length vector indicating desired days of years
%                        of desired download
%  settings            - typical settings structure (see
%                        navsu.internal.initSettings, can also be omitted)
%    .preciseProdDir   - location to put most precise products
%    .galEphCenter     - three letter code indicating the IGS Analysis
%                        Center (AC) to download the precise products from
%    ._____Dir         - lots of other directories that are setup in
%                        initSettings.m given a config file
% OPTIONAL INPUTS:
%  There are various optional inputs depending on the inChoice
%
% OUTPUTS:
%  dayChangei          - days of year when the data has changed locally
%  YearChangei         - years corresponding to dayChangei
% 
% inChoice options:
% 1 :   MGEX Precise Orbit
% 2 :   MGEX Precise Clock
% 3 :   NGA APC
% 4 :   IGS Station Position Solutions
% 5 :   GLONASS Bulletin
% 6 :   Differential code bias
% 7 :   IGS Earth Rotation Parameters
% 8 :   MGEX Observation files
% 9 :   GPS Broadcast Nav Message Files
% 10:   High rate MGEX observation files
% 11:   MGEX Mixed Navigation Files
% 12:   GLONASS Broadcast Nav Message File
% 13:   5 Second CODE GPS clock
% 14:   CNAV data
% 15:   IONEX products
% 16:   BRDM combined multi-GNSS navigation files
% 17:   IGS RINEX 2 Obs Files (non-MGEX core)
% 18:   WAAS NSTB 
% 19:   IGS Real time orbit and clock logs
% 20:   JPL ultra rapid fixing products
% 21:   CODE most recent DCB product
%
% See also: navsu.ftp.ftpFile, navsu.ftp.ftpFileHr, navsu.ftp.websaveFile

if nargin < 4
    % No settings file available
    settings = navsu.internal.initSettings; 
end

cookieFile = settings.cookieFile;
netrcFile = settings.netrcFile;
   
dayChangei  = [];
YearChangei = [];

switch inChoice
    case {1,2} % MGEX precise ephemeris
        if nargin >= 5
            center = varargin{1};
        elseif inChoice == 1
            center = settings.galEphCenter;
        elseif inChoice == 2
            center = settings.galClkCenter;
        end
        
        if inChoice == 1
            exts = {'.sp3','ORB.SP3'};
        else
            exts = {'.clk','CLK.CLK'};
        end
        
        if strcmp(center,'igs') 
            % IGS final is also in a different location
            ftpStruc.destDir      = settings.preciseProdDir;
            ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
            ftpStruc.sourceFormat = '[''/archive/gps/products/'' num2str(gpsWeek) ''/'']';
            ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
            ftpStruc.fileFormat   =  {'[''igs'' num2str(gpsWeek)  num2str(gpsDow) val2 ''.Z'']'};
            ftpStruc.unzipFlag    = 1;
            
            [YearChangei,dayChangei] = navsu.ftp.ftpFile(YearList,dayList,ftpStruc,center,exts{1});
            
        elseif strcmp(center,'igu')
            % IGS final is also in a different location
            ftpStruc.destDir      = settings.preciseProdDir;
            ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
            ftpStruc.sourceFormat = '[''/archive/gps/products/'' num2str(gpsWeek) ''/'']';
            ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
            ftpStruc.fileFormat   =  {'[''igu'' num2str(gpsWeek)  num2str(gpsDow)  ''*''  val2 ''.Z'']'};
            ftpStruc.unzipFlag    = 1;
            
            [YearChangei,dayChangei] = navsu.ftp.ftpFile(YearList,dayList,ftpStruc,center,exts{1});
            
        elseif strcmp(center,'iac')
            
            % GLONASS IAC data is in a different location
            ftpStruc.destDir      = settings.preciseProdDir;
            ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
            ftpStruc.sourceFormat = '[''/pub/glonass/products/'' num2str(gpsWeek) ''/'']';
            ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
            ftpStruc.fileFormat   =  {'[val1 num2str(gpsWeek)  num2str(gpsDow) val2 ''.Z'']'};
            ftpStruc.unzipFlag    = 1;
            
            [YearChangei,dayChangei] = navsu.ftp.ftpFile(YearList,dayList,ftpStruc,center,exts{1});
            
        else
            % RINEX 3 proper naming format
            ftpStruc3.destDir      = settings.preciseProdDir;
            ftpStruc3.ftpSite      = 'https://cddis.nasa.gov';
            ftpStruc3.sourceFormat = '[''/archive/gps/products/mgex/'' num2str(gpsWeek) ''/'']';
            ftpStruc3.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
            ftpStruc3.fileFormat   =  {'[val1 ''0MGXFIN_'' num2str(Year,''%04d'') ''*'' num2str(dayNum,''%03d'') ''*'' val2 ''*'']'};
            ftpStruc3.unzipFlag    = 1;
            
            % Old naming convention
            ftpStruc.destDir      = settings.preciseProdDir;
            ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
            ftpStruc.sourceFormat = '[''/archive/gps/products/mgex/'' num2str(gpsWeek) ''/'']';
            ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
            ftpStruc.fileFormat   =  {'[val1 num2str(gpsWeek)  num2str(gpsDow) val2 ''.Z'']'};
            ftpStruc.unzipFlag    = 1;
            
            center3 = upper(center);
            if strcmp(center3,'COM')
                % they changed their naming code.
                center3 = 'COD';
            end
            
            jdChange = [];
            % try old format
            [YearChange1,dayChange1] = navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile,center,exts{1});
            % try new format!
            [YearChange3,dayChange3] = navsu.ftp.curlFile(YearList,dayList,ftpStruc3,netrcFile,cookieFile,center3,exts{2});
            
            changeOutput = unique(navsu.time.doy2jd([YearChange1; YearChange3],[dayChange1; dayChange3]));
            if ~isempty(changeOutput)
                [dayChangei,YearChangei] = navsu.time.jd2doy(changeOutput);
            end
        end
        
    case 3 % NGA APC Files
        ftpStruc.destDir      = settings.preciseProdDir;
        ftpStruc.ftpSite      = 'ftp.nga.mil';
        ftpStruc.sourceFormat = '[''/pub2/gps/apcpe/'' num2str(Year) ''apc/'']';
        ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''apc'' num2str(gpsWeek) num2str(gpsDow) ''.Z'']' ...
            '[''apc'' num2str(gpsWeek) num2str(gpsDow) ''.exe'']'...
            '[''nga'' num2str(gpsWeek) num2str(gpsDow) ''.apc'']'};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.ftpFile(YearList,dayList,ftpStruc);
        
    case 4 % IGS station position solutions
        % settings.igsStatPos = fslash([baseDir 'Raw Data\IGS_Stations_Final\']);
        ftpStruc.destDir      = settings.preciseProdDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/products/'' num2str(gpsWeek) ''/'']';
        ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
%         ftpStruc.fileFormat   =  {'[''IGS'' num2str(mod(Year,100),''%02i'') ''P'' num2str(woy,''%02i'') ''_all.ssc.Z'']'};
        ftpStruc.fileFormat   =  {'[''IGS'' num2str(mod(Year,100),''%02i'') ''P'' num2str(woy,''%02i'') ''_all.ssc.Z'']' ...
            '[''IGS'' num2str(mod(Year,100),''%02i'') ''P'' num2str(woy,''%02i'') ''.ssc.Z'']'};

        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] = navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
%                 changeOutput = doy2jd(yearChangei,dayChangei);
        
    case 5 % GLONASS bulletin
        ftpStruc.destDir      = settings.gloBulDir;
        ftpStruc.ftpSite      = 'ftp.glonass-iac.ru';
        ftpStruc.sourceFormat = '[''/MCC/BULLETIN/'' num2str(Year) ''/daily/d/'']';
        ftpStruc.destFormat   = '[''/'' int2str(Year) ''/'']';
        ftpStruc.fileFormat   =  {'[''bul_d_'' sprintf(''%02i%02i%02i'',mod(yri,1000),mni,dyi) ''.pdf'']' ;
            '[''bul_d_'' sprintf(''%02i%02i%02i'',mod(yri,1000),mni,dyi) ''.doc'']';};
        ftpStruc.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.ftpFile(YearList,dayList,ftpStruc);
        
    case 6 % IGS differential code bias data
        ftpStruc.destDir      = settings.dcbDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/products/ionex/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.destFormat   = '[int2str(Year) ''\'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''codg'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''i*'']' ;
            '[''casg'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''i*'']' ;};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] = navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
        % MGEX also
        ftpStruc.destDir      = settings.dcbDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/products/bias/'' num2str(Year) ''/'']';
        ftpStruc.destFormat   = '[int2str(Year) ''/'']';
        ftpStruc.fileFormat   =  {'[''DLR0MGXFIN_'' num2str(Year,''%04d'') num2str(max([floor(floor(dayNum/91)*91/10)]),''%02d'') ''*0000_03L_01D_DCB.BSX*'']' };
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] = navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
        % .bia files
        ftpStruc.destDir      = settings.dcbDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/products/mgex/'' num2str(gpsWeek) ''/'']';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''com'' num2str(gpsWeek,''%04d'') num2str(gpsDow) ''.bia.Z'']' };
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
        % RINEX 3 naming convention...
        ftpStruc.destDir      = settings.dcbDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/products/mgex/'' num2str(gpsWeek) ''/'']';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''COD0MGXFIN_'' int2str(Year) num2str(dayNum , ''%03i'') ''0000_01D_01D_OSB.BIA.gz'']' };
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
    case 7 % IGS Earth Rotation Parameters
        ftpStruc.destDir      = settings.erpMgexDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/products/'' num2str(gpsWeek) ''/'']';
        ftpStruc.destFormat   = '[int2str(Year) ''/'']';
        ftpStruc.fileFormat   =  {'[''igs'' num2str(gpsWeek,''%04d'') ''7.erp'']'};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
    case 8 % MGEX Observation Files
        % Optional input of IGS station codes
        if nargin >= 5
            statCodes = varargin{1};
        else
            statCodes = [];
        end
      
        % New RINEX 3 naming convention EOSDIS
        ftpStruc3.destDir      = settings.mgxObsDir;
        ftpStruc3.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc3.sourceFormat = '[''/archive/gps/data/daily/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'' num2str(mod(Year,100),''%02i'') ''d/'' ]';
        ftpStruc3.destFormat   = '[int2str(Year) ''\'' num2str(dayNum,''%03d'') ''/'']';
        
        if isempty(statCodes)
            ftpStruc3.fileFormat   =  {'[''*_30S_MO.crx.gz'']'};
        else
            ftpStruc3.fileFormat = strcat('[''',statCodes,'*_30S_MO.crx.gz'']');
        end
        ftpStruc3.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc3,netrcFile,cookieFile);
        
        % Old RINEX naming convention EOSDIS
        ftpStruc.destDir      = settings.mgxObsDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/data/daily/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'' num2str(mod(Year,100),''%02i'') ''d/'' ]';
        ftpStruc.destFormat   = '[int2str(Year) ''\'' num2str(dayNum,''%03d'') ''/'']';
        
        if isempty(statCodes)
            ftpStruc.fileFormat   =  {'[''*d.Z'']'};
        else
            ftpStruc.fileFormat = strcat('[''',lower(statCodes),'*d.Z'']');
        end
        ftpStruc.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
        
    case 9 % GPS Broadcast Nav Message File
        
        ftpStruc.destDir      = settings.navGpsDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/data/daily/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'' num2str(mod(Year,100),''%02i'') ''n/'' ]';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''*'']'};
        ftpStruc.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
    case 10 % High rate MGEX observation files
        
        error('High rate MGEX observation downloads have not been updated for HTTPS protocols- sorry')
        % Optional input of IGS station codes
        if nargin >= 5
            statCodes = varargin{1};
        else
            statCodes = [];
        end
        
        ftpStruc.destDir      = settings.mgxHrObsDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/data/highrate/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'' num2str(mod(Year,100),''%02i'') ''d/'' num2str(hod,''%02d'') ''/'' ]';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        if isempty(statCodes)
            ftpStruc.fileFormat   =  {'[''*d.Z'']'};
        else
            ftpStruc.fileFormat = strcat('[''',statCodes,'*d.Z'']');
        end
        ftpStruc.unzipFlag    = 0;
        
%         [YearChangei,dayChangei] =navsu.ftp.ftp_file_hr(YearList,dayList,ftpStruc);
        
        ftpStruc3 = ftpStruc;
        if isempty(statCodes)
            ftpStruc3.fileFormat   =  {'[''*o.Z'']'};
        else
            ftpStruc3.fileFormat = strcat('[''',statCodes,'*_01S_MO.crx.gz'']');
        end
        ftpStruc3.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.ftpFileHr(YearList,dayList,ftpStruc3);
        
    case 11 % MGEX Mixed Navigation Files
        ftpStruc.destDir      = settings.navMgxDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/data/daily/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'' num2str(mod(Year,100),''%02i'') ''p/'' ]';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''*MN.rnx*'']'};
        ftpStruc.fileFormat   = {'[''*'']'};
        ftpStruc.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
        % Also download BRDM file
%         ftpHelper(16,YearList,dayList,settings);
        
        
    case 12 % GLONASS Broadcast Nav Message File
        ftpStruc.destDir      = settings.navMgxDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/data/daily/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'' num2str(mod(Year,100),''%02i'') ''g/'' ]';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''*g.Z*'']' '[''*RN.rnx.gz*'']'};
        ftpStruc.fileFormat   =  {'[''*'']'};
        
        ftpStruc.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
    case 13 % CODE 5 Second GPS Clock Data
        ftpStruc.destDir      = settings.preciseProdDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/products/'' num2str(gpsWeek) ''/'']';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        if nargin >= 5
            statCodes = varargin{1};
            ftpStruc.fileFormat   =  {['[''' statCodes ''' num2str(gpsWeek,''%04d'') num2str(gpsDow) ''.clk.Z'']']};
        else
            ftpStruc.fileFormat   =  {'[''cod'' num2str(gpsWeek,''%04d'') num2str(gpsDow) ''.clk_05s.Z'']'};
        end
        
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
        
    case 14 % CNAV data
        ftpStruc.destDir      = settings.navMgxDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/data/campaign/mgex/daily/rinex3/'' num2str(Year) ''/cnav'']';
        ftpStruc.destFormat   = '[int2str(Year) ''/'']';
        ftpStruc.fileFormat   =  {'[''brdx'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''x.Z'']'};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
    case 15 % Iono data
        % IONEX iono map data
        ftpStruc.destDir      = settings.dcbDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/products/ionex/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.destFormat   = '[int2str(Year) ''\'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''codg'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''i*'']' ;
            '[''casg'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''i*'']' ;
            '[''igsg'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''i*'']' ;};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
        % Also pull CODE's klobuchar equivalent
        ftpStruc.destDir      = settings.dcbDir;
        ftpStruc.ftpSite      = 'ftp.aiub.unibe.ch';
        ftpStruc.sourceFormat = '[''/CODE/'' num2str(Year) ''/'' ]';
        ftpStruc.destFormat   = '[int2str(Year) ''\'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''CGIM'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''N.*'']' ;};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.ftpFile(YearList,dayList,ftpStruc);
        
    case 16 % BRDM combined nav files
        ftpStruc.destDir      = settings.navMgxDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gnss/data/campaign/mgex/daily/rinex3/'' num2str(Year) ''/brdm/''  ]';
        ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''brdm'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''p.Z'']'};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
   
        
    case 17 % IGS RINEX 2 Obs Files (non-MGEX core)
        % Optional input of IGS station codes
        if nargin >= 5
            statCodes = varargin{1};
        else
            statCodes = [];
        end
        
        % Old RINEX naming convention
        ftpStruc.destDir      = settings.rinexObsDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/data/daily/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'' num2str(mod(Year,100),''%02i'') ''o/'' ]';
        ftpStruc.destFormat   = '[int2str(Year) ''\'' num2str(dayNum,''%03d'') ''/'']';
        
        if isempty(statCodes)
            ftpStruc.fileFormat   =  {'[''*o.Z'']'};
        else
            ftpStruc.fileFormat = strcat('[''',statCodes,'*o.Z'']');
        end
        ftpStruc.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.ftpFile(YearList,dayList,ftpStruc);
        
    case 18 % WAAS NSTB files
        if nargin >= 5
            statCodes = varargin{1};
        else
            statCodes = [];
        end
        
        if isempty(statCodes)
            ftpStruc.fileFormat   =  {'[''*.gz'']'};
            ftpStruc2017.fileFormat   =  {'[''*.gz'']'};
        else
            ftpStruc.fileFormat = strcat('[''',statCodes,'*.gz'']');
            ftpStruc2017.fileFormat = strcat('[''',statCodes,'*.gz'']');
            
        end
        
        % different years are in different spots, so unfortunately need to
        % separate those
        % transition is november 30, 2017
        inds0 = find((YearList <= 2017 & dayList <= 334) | YearList < 2017);
        inds1 = find(~((YearList <= 2017 & dayList <= 334) | YearList < 2017));

        if ~isempty(inds1)
            ftpStruc.destDir      = settings.waasNstbDir;
            ftpStruc.ftpSite      = 'ftp.nstb.tc.faa.gov';
            ftpStruc.sourceFormat = '[''/pub/NSTB_data/'' num2str(mni,''%02d'') num2str(dyi,''%02d'') num2str(mod(Year,100),''%02d'') ''/''  ]';
            ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
%             ftpStruc.fileFormat   =  {'[''brdm'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''p.Z'']'};
            ftpStruc.unzipFlag    = 1;
            
            [YearChangei,dayChangei] =navsu.ftp.ftpFile(YearList(inds1),dayList(inds1),ftpStruc);
        end
        
        if ~isempty(inds0)
            ftpStruc2017.destDir      = settings.waasNstbDir;
            ftpStruc2017.ftpSite      = 'ftp.nstb.tc.faa.gov';
            ftpStruc2017.sourceFormat = '[''/pub/NSTB_data/000_2017/'' num2str(mni,''%02d'') num2str(dyi,''%02d'') num2str(mod(Year,100),''%02d'') ''/'' ]';
            ftpStruc2017.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
%             ftpStruc2017.fileFormat   =  {'[''brdm'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''p.Z'']'};
            ftpStruc2017.unzipFlag    = 1;
            
            [YearChangei,dayChangei] =navsu.ftp.ftpFile(YearList(inds0),dayList(inds0),ftpStruc2017);
        end
        
    case 19
        % IGS real time products
        ftpStruc.destDir      = settings.preciseProdDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/archive/gps/products/rtpp/'' num2str(gpsWeek) ''/'']';
        ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[''igc'' num2str(gpsWeek)  num2str(gpsDow) ''.sp3.Z'']'};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
    case 20
         % IGS real time products
        ftpStruc.destDir      = settings.preciseProdDir;
        ftpStruc.ftpSite      = 'https://sideshow.jpl.nasa.gov';
        ftpStruc.sourceFormat = '[''/pub/JPL_GNSS_Products/Ultra/'' int2str(Year) ''/'']';
        ftpStruc.destFormat   = '[ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']';
        ftpStruc.fileFormat   =  {'[int2str(Year) ''-'' num2str(mni, ''%02d'') ''-''  num2str(dyi, ''%02d'') ''.pos.gz'']';
            '[int2str(Year) ''-'' num2str(mni, ''%02d'') ''-''  num2str(dyi, ''%02d'') ''.eo.gz'']';
            '[int2str(Year) ''-'' num2str(mni, ''%02d'') ''-''  num2str(dyi, ''%02d'') ''.meta.gz'']';
            '[int2str(Year) ''-'' num2str(mni, ''%02d'') ''-''  num2str(dyi, ''%02d'') ''.pcm.gz'']';
            '[int2str(Year) ''-'' num2str(mni, ''%02d'') ''-''  num2str(dyi, ''%02d'') ''.shadhist.gz'']';
            '[int2str(Year) ''-'' num2str(mni, ''%02d'') ''-''  num2str(dyi, ''%02d'') ''.tdp.gz'']';
            '[int2str(Year) ''-'' num2str(mni, ''%02d'') ''-''  num2str(dyi, ''%02d'') ''.wlpb.gz'']'};
        ftpStruc.unzipFlag    = 1;
        
        [YearChangei,dayChangei] = navsu.ftp.websaveFile(YearList,dayList,ftpStruc);
        
        
    case 21
        
        % 
        ftpStruc.destDir      = settings.dcbDir;
        ftpStruc.ftpSite      = 'https://cddis.nasa.gov';
        ftpStruc.sourceFormat = '[''/gnss/products/bias/'']';
        ftpStruc.destFormat   = '[]';
        ftpStruc.fileFormat   =  {'[''code.bia'']' };
        ftpStruc.unzipFlag    = 0;
        
        [YearChangei,dayChangei] =navsu.ftp.curlFile(YearList,dayList,ftpStruc,netrcFile,cookieFile);
        
    otherwise
        disp('sorry we don''t got that')
end













end