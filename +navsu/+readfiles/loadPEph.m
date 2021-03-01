function [Peph,PFileName,PFileNameFull] = loadPEph(Year, dayNum, settings,FLAG_NO_LOAD,atxData,FLAG_APC_OFFSET,TIME_STRIP)
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
%   .GloPephSource     - this should be set to 'MGEX'- this is old and bad
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

%%
% Adjust in case day number is 0
if dayNum == 0
    Year = Year-1;
    dayNum = navsu.time.yearDays(Year);
end

% Optional flag to not actually load and only pass out filename
if nargin < 4
    FLAG_NO_LOAD = 0;
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
        [Pephi,PFileNamei,PFileNameFulli] = navsu.readfiles.loadPEph(Year(idx), dayNum(idx), settings,FLAG_NO_LOAD,atxData,FLAG_APC_OFFSET,TIME_STRIP);
        
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
    %     [dayNum,Year] = navsu.time.jd2doy(jd);
    gps_day = jd - navsu.time.cal2jd(1980,1,6);
    [yr,mn,dy]=navsu.time.jd2cal(jd);
    [doy,~]=navsu.time.jd2doy(jd);
    
    
    if all(settings.constUse ==  [1 0 0 0 0])
        % NGA APC data
        %precise orbit file
        switch settings.gpsEphCenter
            case 'NGA'
                PfileNameFormat1 = 'NGA%04d%1d.APC';
                PfileNameFormat2 = 'apc%04d%1d';
                PfileNameFormat3 = 'nga%04d%1d.apc';
                
                PpathNameFormat =  [settings.preciseProdDir '/%d/%03d/'];
                
                PFileName = sprintf(PfileNameFormat1, floor((gps_day)/7), mod((gps_day),7));
                if yr >= 2020
                    PFileName = sprintf(PfileNameFormat3, floor((gps_day)/7), mod((gps_day),7));
                    
                elseif yr > 2011
                    PFileName = sprintf(PfileNameFormat2, floor((gps_day)/7), mod((gps_day),7));
                end
                tmp = sprintf(PpathNameFormat, yr,dayNum);
                
                if ~FLAG_NO_LOAD
%                                         Peph = ExpandPeph(navsu.readfiles.readSp3([tmp PFileName]));
                    Peph = (navsu.readfiles.readApc(fullfile([tmp PFileName])));
                    
                    Peph.epochs = ones(Peph.NumSV, 1) * (Peph.GPS_seconds + ...
                        Peph.GPS_week_num*7*24*3600 + ...
                        Peph.Epoch_interval .* (0:Peph.NumEpochs-1));
                    Peph.epochs = Peph.epochs(:);
                end
                PFileNameFull = {[tmp PFileName]};
                PFileName = {PFileName};
                
            case 'IGS'
                PfileNameFormat1 = 'igs%04d%1d.sp3';
                PpathNameFormat =  [settings.preciseProdDir '/%d/%03d/'];
                
                tmp = sprintf(PpathNameFormat, yr,dayNum);
                PFileName = sprintf(PfileNameFormat1, floor((gps_day)/7), mod((gps_day),7));
                
                if ~FLAG_NO_LOAD
                    Peph = ReadSP3([tmp PFileName],0,1,atxData);
                    
                    % If nothing was read, just escape.s
                    if isempty(Peph)
                        return
                    end
                    % Ensure all fields are filled
                    if ~isfield(Peph,'clock_drift')
                        Peph.clock_drift = nan(size(Peph.PRN));
                    end
                    
                    if ~isfield(Peph,'velocity')
                        Peph.velocity = nan(size(Peph.position));
                    end
                    
                    Peph.epochs = ones(Peph.NumSV, 1) * (Peph.GPS_seconds + ...
                        Peph.GPS_week_num*7*24*3600 + ...
                        Peph.Epoch_interval .* (0:Peph.NumEpochs-1));
                    Peph.epochs = Peph.epochs(:);
                end
                PFileNameFull = {[tmp PFileName]};
                PFileName = {PFileName};
                
            case 'jplu'
                % jpl ultra rapids are actually in a different format
                %                 [ int2str(Year) ''/'' num2str(dayNum, ''%03d'') ''/'']
                
                PpathNameFormat =  [settings.preciseProdDir '%d/%03d/'];
                % filename is just year-mn-dy.pos
                PfileNameFormat1 = '%04d-%02d-%02d.pos';
                
                tmp = sprintf(PpathNameFormat, yr,dayNum);
                PFileName = sprintf(PfileNameFormat1, yr, mn,dy);
                
                Peph = navsu.readfiles.readJplPos([tmp PFileName]);
                
            otherwise
                ephCenter = settings.gpsEphCenter;
                
                [Peph,PFileName,PFileNameFull] = loadPephMGEX(ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,1,TIME_STRIP,atxData);
                
        end
        
    elseif all(settings.constUse ==  [0 1 0 0 0])
        GloPephSource = 'MGEX';
        switch GloPephSource
            case 'IGS'
                %precise orbit file
                PfileNameFormat1 = 'igl%04d%1d.sp3';
                PfileNameFormat1 = [settings.gloEphCenter '%04d%1d.sp3'];
                PfileNameFormat2 = 'igl%04d%1d';
                if strcmp(settings.gloEphCenter,'com')
                    % Use CODE MGEX data
                    PpathNameFormat = [settings.rnxMgexPephDir settings.gloEphCenter '/%d/'];
                else
                    PpathNameFormat =  [settings.gloIgsDir '%d/'];
                end
                PFileName = sprintf(PfileNameFormat1, floor((gps_day)/7), mod((gps_day),7));
                tmp = sprintf(PpathNameFormat, yr);
                
                if ~FLAG_NO_LOAD
                    if strcmp(settings.gloEphCenter,'emx') || strcmp(settings.gloEphCenter,'com')
                        Peph = readSp3([tmp PFileName],FLAG_APC_OFFSET,1,2);
                        
                    else
                        Peph = ReadSP3([tmp PFileName],1,1);
                    end
                    % If nothing was read, just escape.s
                    if isempty(Peph)
                        return
                    end
                    % Ensure all fields are filled
                    if ~isfield(Peph,'clock_drift')
                        Peph.clock_drift = nan(size(Peph.PRN));
                    end
                    
                    if ~isfield(Peph,'velocity')
                        Peph.velocity = nan(size(Peph.position));
                    end
                    Peph.epochs = ones(Peph.NumSV, 1) * (Peph.GPS_seconds + ...
                        Peph.GPS_week_num*7*24*3600 + ...
                        Peph.Epoch_interval .* (0:Peph.NumEpochs-1));
                    Peph.epochs = Peph.epochs(:);
                end
                PFileNameFull = {[tmp PFileName]};
                PFileName = {PFileName};
                
            case 'MGEX'
                ephCenter = settings.gloEphCenter;
                
                [Peph,PFileName,PFileNameFull] = loadPephMGEX(ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,2,TIME_STRIP,atxData);
                
        end
        
    elseif all(settings.constUse ==  [0 0 1 0 0])
        % Only have an MGEX option here!
        ephCenter = settings.galEphCenter;
        
        [Peph,PFileName,PFileNameFull] = loadPephMGEX(ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,3,TIME_STRIP,atxData);
        
    elseif all(settings.constUse ==  [0 0 0 1 0])
        % Only have an MGEX option here!
        ephCenter = settings.bdsEphCenter;
        
        [Peph,PFileName,PFileNameFull] = loadPephMGEX(ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,4,TIME_STRIP,atxData) ;
        
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
        
        for cdx = 1:length(settings.constUse)
            if settings.constUse(cdx)
                settings2.constUse = [0 0 0 0 0];
                settings2.constUse(cdx) = 1;
                
                % Call yourself
                [Pephi,PFileNamei,PFileNameFulli] = navsu.readfiles.loadPEph(Year, dayNum, settings2,FLAG_NO_LOAD,atxData,FLAG_APC_OFFSET,TIME_STRIP);
                
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

    function [Peph,PFileName,PFileNameFull] = loadPephMGEX(ephCenter,settings,Year,dayNum,FLAG_NO_LOAD,FLAG_APC_OFFSET,constOut,TIME_STRIP,atxData)
        
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
        PFileName     = [];
        PFileNameFull = [];
        % look for RINEX3 format file first, then fall back on old
        % format if it doesn't exist
        PpathNameFormat =  [settings.preciseProdDir '/%d/%03d/'];
        pathStr = sprintf(PpathNameFormat, Year,dayNum);
        if strcmp(ephCenter,'com')
            center3 = 'COD';
        else
            center3 = upper(ephCenter);
        end
        fname1 = [center3 '0MGXFIN_' num2str(Year,'%04d')  num2str(dayNum,'%03d')];
        % check the directory for a file from that day
        diri = dir(pathStr);
        fileInd = find(~cellfun(@isempty,strfind({diri.name},fname1)) & cellfun(@isempty,strfind({diri.name},'.gz')) & ...
            ~cellfun(@isempty,strfind({diri.name},'ORB.SP3')) );
        
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
            PpathNameFormat =  [settings.preciseProdDir '/%d/%03d/'];
            PFileName = sprintf(PfileNameFormat1, floor((gps_day)/7), mod((gps_day),7));
            tmp = sprintf(PpathNameFormat, yr,dayNum);
        end
        
        if ~FLAG_NO_LOAD
            Peph = navsu.readfiles.readSp3([tmp PFileName],FLAG_APC_OFFSET,1,constOut,atxData);
            
            % Ensure all fields are filled
            if ~isfield(Peph,'clock_drift')
                Peph.clock_drift = nan(size(Peph.PRN));
            end
            
            if ~isfield(Peph,'velocity')
                Peph.velocity = nan(size(Peph.position));
            end
            Peph.epochs = ones(Peph.NumSV, 1) * (Peph.GPS_seconds + ...
                Peph.GPS_week_num*7*24*3600 + ...
                Peph.Epoch_interval .* (0:Peph.NumEpochs-1));
            Peph.epochs = Peph.epochs(:);
            
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
        PFileNameFull = {[tmp PFileName]};
        PFileName = {PFileName};
    end


end


















