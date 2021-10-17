function [Clck, CFileName,CFileNameFull] = loadCFst(Year, dayNum, settings, FLAG_NO_LOAD)
% loadCFst
% DESCRIPTION:
%   Find and parse IGS clock corrections.  The files to be parsed should
%   already exist locally.
%
% INPUTS:
%  Year                - N-length vector of years of desired outputs
%  dayNum              - N-length vector of days of years of desired
%                        outputs
%  settings            - settings structure
%   .constellation     - Desired constellation indicator- can be:
%                        'GPS','GLO','GAL,'BDS','SBAS','MULTI'
%                        If 'MULTI' is chosen, then multiple constellations
%                        will be loaded, and settings.multiConst is
%                        necessary
%   .multiConst        - 1x5 boolean vector indicating desired
%                        constellations, with GPS-GLO-GAL-BDS-SBAS
%                        represented in the respective positions in the
%                        vector
%   .preciseProdDir    - Directory containing precise products- should be
%                        setup in initSettings with a config file
%   .gpsClkCenter      - 3 letter IGS Analysis Center for GPS corrections
%   .gloClkCenter      - 3 letter IGS Analysis Center for GLO corrections
%   .galClkCenter      - 3 letter IGS Analysis Center for GAL corrections
%   .GloPclkSource     - this should be set to 'MGEX'- this is old and bad
%
% OPTIONAL INPUTS:
%  FLAG_NO_LOAD        - True = do not parse the file, just output the name
%                        and locaiton of the local file
%
% OUTPUTS:
%  Clck                - Structure containing parsed precise clock
%                        information
%  CFileName           - Name of precise clock files parsed
%  CFileNameFull       - Name and directory of precise clock files parsed
%
% See also: navsu.ftp.download, navsu.svOrbitClock

% Adjust in case day number is 0
if dayNum == 0
    Year = Year-1;
    dayNum = navsu.time.yearDays(Year);
end

if nargin < 4
    FLAG_NO_LOAD = 0;
end
Clck = [];
if length(dayNum) > 1
    CFileName = {};
    CFileNameFull = {};
    for idx = 1:length(dayNum)
        [Clcki, CFileNamei,filenameFulli] = navsu.readfiles.loadCFst( ...
            Year(idx),dayNum(idx),settings,FLAG_NO_LOAD);
        if ~FLAG_NO_LOAD
            if idx == 1
                Clck = Clcki;
            else
                Clck.Cepochs  = [Clck.Cepochs; Clcki.Cepochs];
                Clck.Cclk     = [Clck.Cclk Clcki.Cclk];
                Clck.Cclk_sig = [Clck.Cclk_sig Clcki.Cclk_sig];

            end
        end
        CFileName = [CFileName CFileNamei];
        CFileNameFull = [CFileNameFull filenameFulli];
    end
    
else
    
    jd = navsu.time.cal2jd(Year,1,0) + dayNum;
    gps_day = jd - navsu.time.cal2jd(1980,1,6);
    [yr,mn,dy]=navsu.time.jd2cal(jd);
    [dayNum,Year]=navsu.time.jd2doy(jd);
    %     dayNum = jd2doy
    
    % Initialize output
    Clck = [];
    
    if all(settings.constUse ==  [1 0 0 0 0])
        switch settings.gpsClkCenter
            case 'IGS'
                %precise clock file
                CfileNameFormat = 'cod%04d%01d.clk_05s';
                CpathNameFormat = '/%d/%03d/';
                
                Clck.Cepochs = [];
                Clck.Cclk = [];
                Clck.Cclk_sig = [];

                for jdx = 1:length(dayNum)
                    % get 30 second clock data from prior and current days
                    CFileName = sprintf(CfileNameFormat, ...
                                        floor(gps_day(jdx)/7), ...
                                        mod(gps_day(jdx), 7));
                    tmp = fullfile(settings.preciseProdDir, ...
                                   sprintf(CpathNameFormat, yr(jdx), dayNum));
                    if ~FLAG_NO_LOAD
                        [Cepochs, Cclk, Cclk_sig] = Read_GPS_05sec_CLK(fullfile(tmp, CFileName),1);
                        Cepochs = Cepochs + 86400*(gps_day(jdx));
                        
                        
                        % Remove known bad points-
                        % Bad clock data on 2015 day 182, PRN 9
                        idx = find(Cepochs > 1119824995 & Cepochs < 1119830395);
                        if ~isempty(idx)
                            Cclk(9,idx) = nan;
                        end
                        
                        Clck.Cepochs  = [ Clck.Cepochs; Cepochs];
                        Clck.Cclk     = [ Clck.Cclk Cclk];
                        Clck.Cclk_sig = [ Clck.Cclk_sig Cclk_sig];
                    end
                end
                
            case 'jplu'
                
                PpathNameFormat = '/%d/%03d/';
                % filename is just year-mn-dy.pos
                PfileNameFormat1 = '%04d-%02d-%02d.tdp';
                
                tmp = fullfile(settings.preciseProdDir, ...
                               sprintf(PpathNameFormat, yr, dayNum));
                CFileName = sprintf(PfileNameFormat1, yr, mn, dy);
                
                % load jpl ultra rapids
                [Cepochs, Cclk, Cclk_sig] = navsu.readfiles.readJplClock( ...
                    fullfile(tmp, CFileName));
               
                Clck.Cepochs  = Cepochs;
                Clck.Cclk     = Cclk;
                Clck.Cclk_sig = Cclk_sig;
               
                Clck.PRNs = (1:size(Clck.Cclk,1))';
                Clck.constInds = 1*ones(size(Clck.PRNs));
            otherwise
                %precise clock file
                clkCenter = settings.gpsClkCenter;
                
                if strcmp(clkCenter,'com')
                    center3 = 'COD';
                else
                    center3 = upper(clkCenter);
                end
                PpathNameFormat = '/%d/%03d/';
                pathStr = fullfile(settings.preciseProdDir, ...
                                   sprintf(PpathNameFormat, yr, dayNum));
                
                fname1 = [center3 '0MGXFIN_' num2str(Year,'%04d')  ...
                          num2str(dayNum,'%03d')];
                % check the directory for a file from that day
                diri = dir(pathStr);
                fileInd = find(contains({diri.name},fname1) ...
                             & ~contains({diri.name},'.gz') ...
                             & contains({diri.name},'.CLK'));
                
                if ~isempty(fileInd)
                    tmp = pathStr;
                    CFileName = diri(fileInd).name;
                else
                    % if the RINEX3 filename is not available, check for the
                    % old one.
                    CfileNameFormat = [clkCenter '%04d%01d.clk'];
                    CpathNameFormat = '/%d/%03d/';
                    CFileName = sprintf(CfileNameFormat, ...
                                        floor((gps_day)/7), ...
                                        mod((gps_day),7));
                    tmp = fullfile(settings.preciseProdDir, ...
                                   sprintf(CpathNameFormat, yr,dayNum));
                end
                
                if ~FLAG_NO_LOAD
                    [Cepochs, Cclk, Cclk_sig] = navsu.readfiles.readRinexClock( ...
                        fullfile(tmp, CFileName), 32, 'G');
                    Cepochs = Cepochs + 86400*(gps_day);
                    
                    %                 if length(Cepochs) == 2880
                    %                    Cepochs = Cepochs(1:10:end);
                    %                    Cclk = Cclk(:,1:10:end);
                    %                    Cclk_sig = Cclk_sig(:,1:10:end);
                    %                 end
                    
                    Clck.Cepochs  = Cepochs;
                    Clck.Cclk     = Cclk;
                    Clck.Cclk_sig = Cclk_sig;

                    
                    Clck.PRNs = (1:size(Clck.Cclk,1))';
                    Clck.constInds = 1*ones(size(Clck.PRNs));
                    
                end
        end
        CFileNameFull = {fullfile(tmp, CFileName)};
        CFileName = {CFileName};
    elseif all(settings.constUse ==  [0 1 0 0 0])
        GloPclkSource = 'MGEX';
        switch GloPclkSource
%             case 'IGS'
%                 %precise clock file
%                 station = 'emx'; % emx/grm
%                 CfileNameFormat = [station '%04d%01d.clk'];
%                 CpathNameFormat =  [settings.preciseProdDir '/%d/%03d/'];
%                 
%                 % get 30 second clock data from prior and current days
%                 CFileName = sprintf(CfileNameFormat, floor((gps_day)/7), mod((gps_day),7));
%                 tmp = sprintf(CpathNameFormat, yr);
%                 if ~FLAG_NO_LOAD
%                     [Cepochs, Cclk, Cclk_sig] = readRinexClock(fullfile(tmp, CFileName),24);
%                     Cepochs = Cepochs + 86400*(gps_day);
%                     
%                     Clck.Cepochs  = Cepochs;
%                     Clck.Cclk     = Cclk;
%                     Clck.Cclk_sig = Cclk_sig;
% 
%                 end
            case 'MGEX'
                clkCenter = settings.gloClkCenter;
                
                PpathNameFormat = '/%d/%03d/';
                pathStr = fullfile(settings.preciseProdDir, ...
                                   sprintf(PpathNameFormat, yr, dayNum));
                
                if strcmp(clkCenter,'com')
                    center3 = 'COD';
                else
                    center3 = upper(clkCenter);
                end
                fname1 = [center3 '0MGXFIN_' num2str(Year,'%04d')  ...
                          num2str(dayNum,'%03d')];
                % check the directory for a file from that day
                diri = dir(pathStr);
                fileInd = find(contains({diri.name},fname1) ...
                            & ~contains({diri.name},'.gz') ...
                            & contains({diri.name},'CLK'));
                
                if ~isempty(fileInd)
                    tmp = pathStr;
                    CFileName = diri(fileInd).name;
                else
                    % if the RINEX3 filename is not available, check for the
                    % old one.
                    CfileNameFormat = [clkCenter '%04d%01d.clk'];
                    CpathNameFormat = '/%d/%03d/';
                    CFileName = sprintf(CfileNameFormat, ...
                                        floor((gps_day)/7), ...
                                        mod((gps_day),7));
                    tmp = fullfile(settings.preciseProdDir, ...
                                   sprintf(CpathNameFormat, yr, dayNum));
                end
                
                if ~FLAG_NO_LOAD
                    [Cepochs, Cclk, Cclk_sig] = navsu.readfiles.readRinexClock( ...
                        fullfile(tmp, CFileName), 24, 'R');
                    Cepochs = Cepochs + 86400*(gps_day);
                    
                    %                 if length(Cepochs) == 2880
                    %                    Cepochs = Cepochs(1:10:end);
                    %                    Cclk = Cclk(:,1:10:end);
                    %                    Cclk_sig = Cclk_sig(:,1:10:end);
                    %                 end
                    
                    Clck.Cepochs  = Cepochs;
                    Clck.Cclk     = Cclk;
                    Clck.Cclk_sig = Cclk_sig;

                    
                    Clck.PRNs = (1:size(Clck.Cclk,1))';
                    Clck.constInds = 2*ones(size(Clck.PRNs));
                    
                end
        end
        CFileNameFull = {fullfile(tmp, CFileName)};
        CFileName = {CFileName};
    elseif all(settings.constUse ==  [0 0 1 0 0])
        clkCenter = settings.galClkCenter;
        
        PpathNameFormat = '/%d/%03d/';
        pathStr = fullfile(settings.preciseProdDir, ...
                           sprintf(PpathNameFormat, yr, dayNum));
        
        if strcmp(clkCenter,'com')
            center3 = 'COD';
        else
            center3 = upper(clkCenter);
        end
        fname1 = [center3 '0MGXFIN_' num2str(Year,'%04d')  num2str(dayNum,'%03d')];
        % check the directory for a file from that day
        diri = dir(pathStr);
        fileInd = find(contains({diri.name}, fname1) ...
                    & ~contains({diri.name}, '.gz') ...
                    & contains({diri.name}, 'CLK')   );
        
        if ~isempty(fileInd)
            tmp = pathStr;
            CFileName = diri(fileInd).name;
        else
            % if the RINEX3 filename is not available, check for the
            % old one.
            CfileNameFormat = [clkCenter '%04d%01d.clk'];
            CpathNameFormat = '/%d/%03d/';
            CFileName = sprintf(CfileNameFormat, floor((gps_day)/7), mod((gps_day),7));
            tmp = fullfile(settings.preciseProdDir, ...
                           sprintf(CpathNameFormat, yr,dayNum));
        end
        
        if ~FLAG_NO_LOAD
            [Cepochs, Cclk, Cclk_sig] = navsu.readfiles.readRinexClock(fullfile(tmp, CFileName),36,'E');
            Cepochs = Cepochs + 86400*(gps_day);
            
            Clck.Cepochs  = Cepochs;
            Clck.Cclk     = Cclk;
            Clck.Cclk_sig = Cclk_sig;

            
            Clck.PRNs = (1:size(Clck.Cclk,1))';
            Clck.constInds = 3*ones(size(Clck.PRNs));
        end
        CFileNameFull = {fullfile(tmp, CFileName)};
        CFileName = {CFileName};
        
    elseif all(settings.constUse ==  [0 0 0 1 0])
        clkCenter = settings.bdsClkCenter;
        
        PpathNameFormat = '/%d/%03d/';
        pathStr = fullfile(settings.preciseProdDir, ...
                           sprintf(PpathNameFormat, yr, dayNum));
                       
        if strcmp(clkCenter,'com')
            center3 = 'COD';
        else
            center3 = upper(clkCenter);
        end
        fname1 = [center3 '0MGXFIN_' num2str(Year,'%04d')  num2str(dayNum,'%03d')];
        % check the directory for a file from that day
        diri = dir(pathStr);
        fileInd = find(contains({diri.name}, fname1) ...
                     & ~contains({diri.name}, '.gz') ...
                     & ~contains({diri.name}, 'CLK'));
        
        if ~isempty(fileInd)
            tmp = pathStr;
            CFileName = diri(fileInd).name;
        else
            % if the RINEX3 filename is not available, check for the
            % old one.
            CfileNameFormat = [clkCenter '%04d%01d.clk'];
            CpathNameFormat = '/%d/%03d/';
            CFileName = sprintf(CfileNameFormat, floor((gps_day)/7), mod((gps_day),7));
            tmp = fullfile(settings.preciseProdDir, ...
                           sprintf(CpathNameFormat, yr, dayNum));
        end
        
        if ~FLAG_NO_LOAD
            [Cepochs, Cclk, Cclk_sig] =  navsu.readfiles.readRinexClock(fullfile(tmp, CFileName),35,'C');
            Cepochs = Cepochs + 86400*(gps_day);
            
            Clck.Cepochs  = Cepochs;
            Clck.Cclk     = Cclk;
            Clck.Cclk_sig = Cclk_sig;
            
            Clck.PRNs = (1:size(Clck.Cclk,1))';
            Clck.constInds = 4*ones(size(Clck.PRNs));
        end
        CFileNameFull = {fullfile(tmp, CFileName)};
        CFileName = {CFileName};
        
    elseif all(settings.constUse ==  [0 0 0 0 1])
        % SBAS case currently not covered. Just return empty outputs to
        % avoid infinite recursion.
        
        Clck.Cepochs = [];
        Clck.Cclk = [];
        Clck.Cclk_sig = [];
        
        CFileName = '';
        CFileNameFull = '';
        
    else
        settings2 = settings;
        
        CFileName = {};
        
        PRNs = [];
        constInds = [];
        CFileNameFull = {};
        
        for cdx = 1:length(settings.constUse)
            if settings.constUse(cdx)
                settings2.constUse = [0 0 0 0 0];
                settings2.constUse(cdx) = 1;
                
                % Call yourself
                [Cdatai,CFileNamei,CFileNameFulli] = navsu.readfiles.loadCFst( ...
                    Year, dayNum, settings2,FLAG_NO_LOAD);
                
                if ~FLAG_NO_LOAD
                    if isempty(PRNs)
                        PRNs = (1:size(Cdatai.Cclk,1))';
                        constInds = cdx*ones(size(PRNs));
                        
                        Clck.Cclk = Cdatai.Cclk;
                        Clck.Cepochs = Cdatai.Cepochs;
                        Clck.Cclk_sig = Cdatai.Cclk_sig;

                    else
                        PRNs = [PRNs; (1:size(Cdatai.Cclk,1))'];
                        constInds = [constInds; cdx*ones(size(Cdatai.Cclk,1),1)];
                        
                        Clck.Cclk = [Clck.Cclk; Cdatai.Cclk];
                        Clck.Cclk_sig = [Clck.Cclk_sig; Cdatai.Cclk_sig];
                    end
                end
                
                CFileName = [CFileName CFileNamei];
                CFileNameFull = [CFileNameFull CFileNameFulli];
            end
        end
        Clck.PRNs = PRNs;
        Clck.constInds = constInds;
        
    end
    
end
end




