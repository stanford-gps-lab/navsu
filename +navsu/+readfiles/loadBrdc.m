function [eph,BFileName,BFileNameFull] = loadBrdc(Year, dayNum, settings,varargin)
%% Load broadcast navigation message

% Adjust in case day number is 0
if dayNum == 0
    Year = Year-1;
    dayNum = YearDays(Year);
end

p = inputParser;
p.addParameter('outFormat', 'struct'); % 'struct' or 'array'
p.addParameter('FLAG_NO_LOAD',0);


% parse the results
parse(p, varargin{:});
res = p.Results;
outFormat      = res.outFormat;
FLAG_NO_LOAD   = res.FLAG_NO_LOAD;

BpathName = [];
BFileName = [];

eph = struct('gps',[],'glo',[],'gal',[],'bds',[],'qzss',[],'iono',[],'leapSecond',[]);

if length(dayNum) > 1
    for idx = 1:length(dayNum)
        [ephi,BFileNamei,BFileNameFulli] = navsu.readfiles.loadBrdc(Year(idx), dayNum(idx), ...
            settings,'FLAG_NO_LOAD',FLAG_NO_LOAD,'outFormat','array');
        
        % Add everything to the existing array
        eph.gps = cat(2,eph.gps,ephi.gps);
        eph.glo = cat(2,eph.glo,ephi.glo);
        eph.gal = cat(2,eph.gal,ephi.gal);
        eph.bds = cat(2,eph.bds,ephi.bds);        
        eph.qzss = cat(2,eph.qzss,ephi.qzss);
        eph.iono = ephi.iono;
        eph.leapSecond = ephi.leapSecond;

        BFileName{idx} = BFileNamei;
        BFileNameFull{idx} = BFileNameFulli;
    end
    
else
    % single day- can be multi-constellation!
    jd = navsu.time.cal2jd(Year,1,0) + dayNum;
    gps_day = jd - navsu.time.cal2jd(1980,1,6);
    [yr,mn,dy]= navsu.time.jd2cal(jd);
    [doy,~]= navsu.time.jd2doy(jd);
    
    % Default multi-GNSS file
    BfileNameFormat = '%4d/%03d/brdm%03d0.%02dp';
    brdmPath = settings.navMgxDir;
    brdmFilename = sprintf(BfileNameFormat, yr, doy, doy, mod(yr, 100));
    
    if settings.constUse(1)
        % GPS
        if isfield(settings,'gpsNavSource') && strcmp(settings.gpsNavSource,'s')
            BfileNameFormat = '%4d/sugl%03d0.%02dn';
            BpathName = settings.suglGpsDir;
            BFileName = sprintf(BfileNameFormat, yr, doy, mod(yr, 100));
        else
            % MGEX file
            BFileName = brdmFilename;
            BpathName = brdmPath;
        end
        constellations = navsu.readfiles.initConstellation(1,0,0,0,0);
        
        BFileNameFull{1} = fullfile([BpathName BFileName]);
        if ~FLAG_NO_LOAD
            ephi = navsu.readfiles.loadRinexNav(fullfile([BpathName BFileName]),'constellations',constellations, 'outFormat','array');
            eph.gps = ephi.gps;
        end
    end
    
    if settings.constUse(3)
        % Galileo
        if strcmp(settings.galNavSource,'c') % BKG product
            BfileNameFormat = '%4d/%03d/BRDC00WRD_R_%03d%03d0000_01D_MN.rnx';
            BpathName = settings.rnxMgexNavDir;
            BFileName = sprintf(BfileNameFormat, yr, doy, yr,doy);
        elseif strcmp(settings.galNavSource ,'m')  % default to the TUM/DLR files
            BfileNameFormat = '%4d/%03d/brdm%03d0.%02dp';
            BpathName = settings.rnxMgexNavDir;
            BFileName = sprintf(BfileNameFormat, yr, doy, doy, mod(yr, 100));
        elseif strcmp(settings.galNavSource, 's') % SUGL file lol
            BfileNameFormat = '%4d/sugl%03d0.%02dl';
            BpathName = settings.suglGalDir;
            
            BFileName = sprintf(BfileNameFormat, yr, doy, mod(yr, 100));
        end
        BFileNameFull{1} = fullfile([BpathName BFileName]);

        constellations = navsu.readfiles.initConstellation(0, 0, 1, 0, 0, 0);
        if ~FLAG_NO_LOAD
            ephi = navsu.readfiles.loadRinexNav(fullfile([BpathName BFileName]),'constellations',constellations, 'outFormat','array');
            eph.gal = ephi.gal;
        end
    end
    
        
end

if strcmp(outFormat,'struct')
    eph.gps  = navsu.readfiles.ephArray2Struct(eph.gps',[],eph.leapSecond,'GPS');
    eph.glo  = navsu.readfiles.ephArray2StructGlonass(eph.glo',[],eph.leapSecond);
    eph.gal  = navsu.readfiles.ephArray2Struct(eph.gal',[],eph.leapSecond,'GAL');
    eph.bds  = navsu.readfiles.ephArray2Struct(eph.bds',[],eph.leapSecond,'BDS');
    eph.qzss = navsu.readfiles.ephArray2Struct(eph.qzss',[],eph.leapSecond,'QZSS');
    eph.iono = eph.iono;
end

end

