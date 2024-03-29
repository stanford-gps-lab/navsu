function [Peph,PFileName,PFileNameFull] = loadPEph(Year, dayNum, settings, ...
    FLAG_NO_LOAD, atxData, FLAG_APC_OFFSET, TIME_STRIP)
% loadPEph
% DESCRIPTION:
% Find and parse IGS clock corrections.  The files to be parsed should
% already exist locally.
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
%   .gpsEphCenter      - 3 letter IGS Analysis Center for GPS corrections
%                        This can be 'NGA','IGS', or one of the MGEX AC's.
%                        'IGS' indicates the IGS final solution
%   .gloEphCenter      - 3 letter IGS Analysis Center for GLO corrections
%   .galEphCenter      - 3 letter IGS Analysis Center for GAL corrections
%
% OPTIONAL INPUTS:
%  FLAG_NO_LOAD        - True = do not parse the file, just output the name
%                        and location of the local file
%  atxData             - Structure of IGS antenna phase center information
%  FLAG_APC_OFFSET     - True = displace the interpolated positions by their
%                        antenna phase center offset. Default = true.
%
% OUTPUTS:
%  Peph                - Structure containing parsed precise orbit
%                        information
%  PFileName           - Name of precise orbit files parsed
%  PFileNameFull       - Name and directory of precise orbit files parsed
%
% See also: navsu.ftp.download, navsu.svOrbitClock

if isscalar(Year)
    Year = repmat(Year, size(dayNum));
end

%%
total_days = navsu.time.YearDays(Year);
% Adjust in case day number is 0
day0 = dayNum == 0;
if any(day0)
    Year(day0) = Year(day0)-1;
    dayNum(day0) = total_days(day0);
end

% Adjust in case day number is 366 or 367 for year start/end
day366 = dayNum == total_days + 1;
if any(day366)
    Year(day366) = Year(day366) + 1;
    dayNum(day366) = 1;
end

% Optional flag to not actually load and only pass out filename
if nargin < 4
    FLAG_NO_LOAD = false;
end

% optional antenna phase center structure
if nargin < 5
    atxData = [];
end

if nargin < 6
    FLAG_APC_OFFSET = 1;
end
if nargin < 7
    TIME_STRIP = 0;
end

if length(dayNum) > 1
    PFileName = {}; PFileNameFull = {};
    for idx = 1:length(dayNum)
        [Pephi,PFileNamei,PFileNameFulli] = navsu.readfiles.loadPEph( ...
            Year(idx), dayNum(idx), settings,FLAG_NO_LOAD,atxData,FLAG_APC_OFFSET,TIME_STRIP);
        
        if idx == 1
            Peph = Pephi;
        else
            Peph.PRN           = [Peph.PRN; Pephi.PRN];
            Peph.clock_bias    = [Peph.clock_bias; Pephi.clock_bias];
            Peph.position      = [Peph.position; Pephi.position];
            Peph.Event         = [Peph.Event; Pephi.Event];
            Peph.clock_drift   = [Peph.clock_drift; Pephi.clock_drift];
            Peph.velocity      = [Peph.velocity; Pephi.velocity];
            Peph.epochs        = [Peph.epochs; Pephi.epochs];
            if isfield(Peph,'constellation')
                Peph.constellation = [Peph.constellation; Pephi.constellation];
            end
        end
        PFileName = [PFileName PFileNamei];
        PFileNameFull = [PFileNameFull PFileNameFulli];
    end
else
    Peph = [];
    
    jd = navsu.time.cal2jd(Year,1,0) + floor(dayNum);
    % adjust Year and dayNum in case of rollover
    gps_day = jd - navsu.time.cal2jd(1980,1,6);
    [yr,mn,dy] = navsu.time.jd2cal(jd);
    
    
    if all(settings.constUse ==  [1 0 0 0 0])
        
        if any(strcmp(settings.gpsEphCenter, {'NGA', 'IGS', 'jplu'}))
            % individual processing to get PFileName

            switch settings.gpsEphCenter
                case 'NGA'

                    if yr >= 2020
                        PfileNameFormat = 'nga%04d%1d.apc';
                    elseif yr > 2011
                        PfileNameFormat = 'apc%04d%1d';
                    else
                        PfileNameFormat = 'NGA%04d%1d.APC';
                    end
                    PFileName = sprintf(PfileNameFormat, ...
                                        floor(gps_day/7), ...
                                        mod(gps_day, 7));

                case 'IGS'
                    PfileNameFormat = 'igs%04d%1d.sp3';

                    PFileName = sprintf(PfileNameFormat, ...
                                        floor(gps_day/7), ...
                                        mod(gps_day, 7));

                case 'jplu'
                    % jpl ultra rapids are actually in a different format
                    %                 [ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']

                    % filename is just year-mn-dy.pos
                    PfileNameFormat1 = '%04d-%02d-%02d.pos';

                    PFileName = sprintf(PfileNameFormat1, yr, mn, dy);

            end

            % now compile file name:
            PpathNameFormat = '/%d/%03d/'; % can be set for all

            tmp = sprintf(PpathNameFormat, yr, dayNum);
            fullFilePathName = fullfile(settings.preciseProdDir, tmp, PFileName);

            PFileName = {PFileName};
            PFileNameFull = {fullFilePathName};
        end
        
        % now actually download
        switch settings.gpsEphCenter
            case 'NGA'
                if ~FLAG_NO_LOAD
%                     Peph = navsu.readfiles.readSp3(fullFilePathName);
                    Peph = navsu.readfiles.readApc(fullFilePathName);
                    
                end
                
            case 'IGS'
                
                if ~FLAG_NO_LOAD
                    Peph = navsu.readfiles.readSp3(fullFilePathName, 0, 1);
                    
                    % If nothing was read, just escape.s
                    if isempty(Peph)
                        return
                    end
                    % Ensure all fields are filled
                    if ~isfield(Peph,'clock_drift')
                        Peph.clock_drift = nan(size(Peph.PRN));
                    end
                    
                end
                
            case 'jplu'
                
                Peph = navsu.readfiles.readJplPos(fullFilePathName);
                
            otherwise
                ephCenter = settings.gpsEphCenter;
                
                % this function does it all
                [Peph,PFileName,PFileNameFull] = loadPephMGEX( ...
                    ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,1,TIME_STRIP,atxData);
                
        end
        
    elseif all(settings.constUse == [0 1 0 0 0])
        GloPephSource = 'MGEX';
        switch GloPephSource
            case 'IGS'
                %precise orbit file
                PfileNameFormat1 = 'igl%04d%1d.sp3';
                PfileNameFormat1 = [settings.gloEphCenter '%04d%1d.sp3'];
                PFileName = sprintf(PfileNameFormat1, floor((gps_day)/7), mod((gps_day),7));
                
                if strcmp(settings.gloEphCenter,'com')
                    % Use CODE MGEX data
                    PpathNameFormat = [settings.rnxMgexPephDir settings.gloEphCenter '/%d/'];
                else
                    PpathNameFormat =  [settings.gloIgsDir '%d/'];
                end
                tmp = sprintf(PpathNameFormat, yr);
                
                fullFilePathName = fullfile(tmp, PFileName);

                if ~FLAG_NO_LOAD
                    if strcmp(settings.gloEphCenter,'emx') || strcmp(settings.gloEphCenter,'com')
                        Peph = navsu.readfiles.readSp3(fullFilePathName,FLAG_APC_OFFSET,1,2);
                        
                    else
                        Peph = navsu.readfiles.readSp3(fullFilePathName,1,1);
                    end
                    % If nothing was read, just escape.s
                    if isempty(Peph)
                        return
                    end
                    % Ensure all fields are filled
                    if ~isfield(Peph,'clock_drift')
                        Peph.clock_drift = nan(size(Peph.PRN));
                    end
                    
                end
                PFileNameFull = {fullFilePathName};
                PFileName = {PFileName};
                
            case 'MGEX'
                ephCenter = settings.gloEphCenter;
                
                [Peph,PFileName,PFileNameFull] = loadPephMGEX( ...
                    ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,2,TIME_STRIP,atxData);
                
        end
        
    elseif all(settings.constUse == [0 0 1 0 0])
        % Only have an MGEX option here!
        ephCenter = settings.galEphCenter;
        
        [Peph,PFileName,PFileNameFull] = loadPephMGEX( ...
            ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,3,TIME_STRIP,atxData);
        
    elseif all(settings.constUse == [0 0 0 1 0])
        % Only have an MGEX option here!
        ephCenter = settings.bdsEphCenter;
        
        [Peph,PFileName,PFileNameFull] = loadPephMGEX( ...
            ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,4,TIME_STRIP,atxData) ;
    elseif all(settings.constUse == [0 0 0 0 1])
        % there is a mismatch here! svOrbitClock considers this QZSS. Other
        % functions consider this SBAS and QZSS is #6. Either way 
        % navsu.readfiles.readSp3 can not handle it. But this case
        % differentiation is necessary to avoid an infinite recursion.
        warning('Currently unable to parse ephemeris for constellation #5.')
        Peph.PRN           = [];
        Peph.clock_bias    = [];
        Peph.clock_drift   = [];
        Peph.position      = [];
        Peph.velocity      = [];
        Peph.Event         = [];
        Peph.epochs        = [];
        Peph.constellation = [];
        PFileName = '';
        PFileNameFull = '';
    else
        % Multi-GNSS
        PRN           = [];
        clock_bias    = [];
        clock_drift   = [];
        position      = [];
        velocity      = [];
        Event         = [];
        epochs        = [];
        constellation = [];
        PFileName = {};
        PFileNameFull = {};
        
        settings2 = settings;
        
        for cdx = find(settings.constUse)
            settings2.constUse = [0 0 0 0 0];
            settings2.constUse(cdx) = 1;

            % Call yourself
            [Pephi,PFileNamei,PFileNameFulli] = navsu.readfiles.loadPEph( ...
                Year, dayNum, settings2,FLAG_NO_LOAD,atxData,FLAG_APC_OFFSET,TIME_STRIP);

            if ~FLAG_NO_LOAD
                PRN           = [PRN; Pephi.PRN];
                clock_bias    = [clock_bias; Pephi.clock_bias];
                clock_drift   = [clock_drift; Pephi.clock_drift];
                position      = [position; Pephi.position];
                velocity      = [velocity; Pephi.velocity];
                Event         = [Event; Pephi.Event];
                epochs        = [epochs; Pephi.epochs];
                constellation = [constellation; cdx*ones(size(Pephi.epochs))];
            end
            PFileName      = [PFileName PFileNamei];
            PFileNameFull = [PFileNameFull PFileNameFulli];
        end
        
        Peph.PRN           = PRN;
        Peph.clock_bias    = clock_bias;
        Peph.clock_drift   = clock_drift;
        Peph.position      = position;
        Peph.velocity      = velocity;
        Peph.Event         = Event;
        Peph.epochs        = epochs;
        Peph.constellation = constellation;
    end
end

    function [Peph,PFileName,PFileNameFull] = loadPephMGEX( ...
            ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,constOut,TIME_STRIP,atxData)
        
        % day of year can sometimes be fractional for ultra rapids- need to
        % clean it up
        dayFrac = rem(dayNum,1);
        dayNum = floor(dayNum);
        
        Peph.PRN           = [];
        Peph.clock_bias    = [];
        Peph.clock_drift   = [];
        Peph.position      = [];
        Peph.velocity      = [];
        Peph.Event         = [];
        Peph.epochs        = [];
        Peph.constellation = [];
        % look for RINEX3 format file first, then fall back on old
        % format if it doesn't exist
        PpathNameFormat = '/%d/%03d/';
        pathStr = fullfile(settings.preciseProdDir, sprintf(PpathNameFormat, Year,dayNum));
        if strcmp(ephCenter,'com')
            center3 = 'COD';
        else
            center3 = upper(ephCenter);
        end
        fname1 = [center3 '0MGXFIN_' num2str(Year,'%04d')  num2str(dayNum,'%03d')];
        % check the directory for a file from that day
        diri = dir(pathStr);
        fileInd = find(contains({diri.name}, fname1) ...
                     & ~contains({diri.name}, '.gz') ...
                     & contains({diri.name}, 'ORB.SP3'));
        
        if ~isempty(fileInd)
            tmp = pathStr;
            PFileName = diri(fileInd).name;
        else
            % if the RINEX3 filename is not available, check for the
            % old one.
            
            if strcmp(ephCenter,'igu')
                fileHour = floor(dayFrac/0.25)*6;
                PfileNameFormat1 = [ephCenter '%04d%1d_' num2str(fileHour,'%02i') '.sp3'];
                
            else
                PfileNameFormat1 = [ephCenter '%04d%1d.sp3'];
            end
            PpathNameFormat = '/%d/%03d/';
            PFileName = sprintf(PfileNameFormat1, floor((gps_day)/7), mod((gps_day),7));
            tmp = fullfile(settings.preciseProdDir, sprintf(PpathNameFormat, yr,dayNum));
        end
        
        fullFileName = fullfile(tmp, PFileName);
        
        if ~FLAG_NO_LOAD
            Peph = navsu.readfiles.readSp3(fullFileName,FLAG_APC_OFFSET,1,constOut,atxData);
            
            % Ensure all fields are filled
            if ~isfield(Peph,'clock_drift')
                Peph.clock_drift = nan(size(Peph.PRN));
            end
            
            
            Peph.Event(isnan(Peph.Event)) = 1;
            
            
            % strip off epochs that aren't on correct day
            if TIME_STRIP
                doysi = floor(navsu.time.jd2doy(navsu.time.epochs2jd(Peph.epochs)));
                indsRemove = find(doysi ~= dayNum);
                if ~isempty(indsRemove)
                    Peph.PRN(indsRemove) = [];
                    Peph.clock_bias(indsRemove) = [];
                    Peph.position(indsRemove,:) = [];
                    Peph.Event(indsRemove) = [];
                    Peph.clock_drift(indsRemove) = [];
                    Peph.velocity(indsRemove,:) = [];
                    Peph.epochs(indsRemove) = [];
                end
            end
            
            Peph.constellation = constOut*ones(size(Peph.epochs));
        end
        PFileNameFull = {fullFileName};
        PFileName = {PFileName};
    end


end
