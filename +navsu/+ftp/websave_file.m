function  [YearChange,dayChange] = websave_file(YearList,dayList,ftpStruc)

ftpSite      = ['' ftpStruc.ftpSite];
sourceFormat = ftpStruc.sourceFormat;
destDir      = ftpStruc.destDir;
destFormat   = ftpStruc.destFormat;
fileFormat   = ftpStruc.fileFormat;
unzipFlag    = ftpStruc.unzipFlag;

YearChange = [];
dayChange  = [];

for ddx = 1:length(YearList)
    dayNum = dayList(ddx);
    Year = YearList(ddx);
    
    jdi = navsu.time.doy2jd(Year,dayNum);
    [yri,mni,dyi] = navsu.time.jd2cal(jdi);
    [gpsWeek,gpsTow] = navsu.time.jd2gps(jdi);
    gpsDow = floor(gpsTow/86400);
    
    % Initial week of year
    [gpsWeek0,tow0] = navsu.time.jd2gps(navsu.time.cal2jd(Year,1,1));
    woy = gpsWeek-gpsWeek0+1;
    
    target_dir = [destDir eval(destFormat)];
    
    if ~exist(target_dir,'dir')
       mkdir(target_dir) 
    end
    
    sourceDir = [ftpSite eval(sourceFormat)];
    
    for fdx = 1:length(ftpStruc.fileFormat)
        % Download the file for each day
        desiredName = eval(fileFormat{fdx});
        
        websave([target_dir desiredName],[sourceDir desiredName]);
        
        if unzipFlag
            navsu.readfiles.unzipFile([target_dir desiredName]);
        end
        
    end
end














end