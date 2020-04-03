function [dcbData, filesRead] = loadDcb(YearList,dayList,settings,FLAG_NO_LOAD,SOURCE_SELECT,statCode)

if nargin < 4
    FLAG_NO_LOAD = 0;
end

if nargin < 5
    SOURCE_SELECT = 4;
end

if nargin < 6
    statCode = [];
end

dayStartEpochs = utility.time.jd2epochs(utility.time.doy2jd(YearList,dayList));

dcbData = [];
filesRead = {};
for ddx = 1:length(YearList)
    
    Year = YearList(ddx);
    dayNum = dayList(ddx);
    [gpsWeek,gpsTow] = utility.time.jd2gps(utility.time.doy2jd(Year,dayNum));
    gpsDow = floor(gpsTow/86400);
    
    
    switch SOURCE_SELECT
        case 1
            % this is not really to be used right now
            
            % DCB File format
            DfileNameFormat = '%4d/%03d/CAS0MGXRAP_%4d%03d0000_01D_01D_DCB.BSX';
            
            % File location
            DpathName = settings.dcbMgexDir;
            
            % Filename
            DFileName = sprintf(DfileNameFormat, Year, dayNum,Year, dayNum);
            
            % Load the file
            dcbDatai = parseDcbBsxFile([DpathName DFileName]);
        case 2
            if dayNum >= 364
                dayNum = 363;
            end
            
            % DLR product
            ftpStruc.destDir      = settings.dcbDir;
            ftpStruc.ftpSite      = 'cddis.gsfc.nasa.gov';
            ftpStruc.sourceFormat = '[''/gnss/products/bias/'' num2str(Year) ''/'']';
            ftpStruc.destFormat   = '[int2str(Year) ''/'']';
            ftpStruc.fileFormat   =  {'[''DLR0MGXFIN_'' num2str(Year,''%04d'') num2str(floor(floor(dayNum/91)*91/10),''%02d'') ''*0000_03L_01D_DCB.BSX*'']' };
            filenamei = [settings.dcbDir eval(ftpStruc.destFormat) eval(ftpStruc.fileFormat{1})];
            
            filenamei = filenamei(1:end-1);
            %         filenamei(strfind(filenamei,'*')) = [];
            
            [fold,filei,exti] = fileparts(filenamei);
            
            files = dir(fileparts(filenamei));
            found = 0;
            for idx = 1:length(files)
                desiredName = regexptranslate('wildcard',[filei exti]);
                
                isthisit = ~isempty(regexp(files(idx).name,desiredName,'once')) && isempty(strfind(files(idx).name,'.gz')) ;
                
                if isthisit
                    filenamei = [fold '/' files(idx).name];
                    found = 1;
                    break
                end
            end
            
            dcbDatai = [];
            if found && isempty(utility.strFindCell(filesRead,filenamei))
                filesRead = [filesRead {filenamei}];
                if ~FLAG_NO_LOAD
                    dcbDatai = utility.readfiles.parseDcbBsxFile(filenamei);
                    
                    % Only include the current day
                    indsInclude = find(ismember(dcbDatai.startEpoch,dayStartEpochs));
                    
                    fields = fieldnames(dcbDatai);
                    for idx = 1:length(fields)
                        if length(dcbDatai.(fields{idx})) > 5
                            dcbDatai.(fields{idx}) = dcbDatai.(fields{idx})(indsInclude);
                        end
                    end
                end
            end
            
        case 3
            % this is not really to be used right now
            % nor this
            ftpStruc.destDir      = settings.dcbDir;
            ftpStruc.ftpSite      = 'cddis.gsfc.nasa.gov';
            ftpStruc.sourceFormat = '[''/gnss/products/ionex/'' num2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
            ftpStruc.destFormat   = '[int2str(Year) ''\'' num2str(dayNum,''%03d'') ''/'']';
            ftpStruc.fileFormat   =  {'[''codg'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''i*'']' ;
                '[''casg'' num2str(dayNum,''%03d'') ''0.'' num2str(mod(Year,100),''%02d'') ''i*'']' ;};
            ftpStruc.unzipFlag    = 1;
            
            
            filenamei = [settings.dcbDir eval(ftpStruc.destFormat) eval(ftpStruc.fileFormat{2})];
            
            filenamei(strfind(filenamei,'*')) = [];
            
            dcbDatai = utility.readfiles.readIONEX(filenamei);
        case 4
            % Daily CODE DCB estimates.
            ftpStruc.destDir      = settings.dcbDir;
            ftpStruc.ftpSite      = 'cddis.gsfc.nasa.gov';
            ftpStruc.sourceFormat = '[''/pub/gps/products/mgex/'' num2str(gpsWeek) ''/'']';
            ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
            ftpStruc.fileFormat   =  {'[''com'' num2str(gpsWeek,''%04d'') num2str(gpsDow) ''.bia.Z'']' };
            ftpStruc.unzipFlag    = 1;
            filenamei = [settings.dcbDir eval(ftpStruc.destFormat) eval(ftpStruc.fileFormat{1})];
            filenamei = filenamei(1:end-2);
            
            if exist(filenamei,'file')
                dcbDatai = parseDcbBsxFile(filenamei,0);
            else
                % Need to check for RINEX3 naming convention!
                ftpStruc.destDir      = settings.dcbDir;
                ftpStruc.ftpSite      = 'cddis.gsfc.nasa.gov';
                ftpStruc.sourceFormat = '[''/pub/gps/products/mgex/'' num2str(gpsWeek) ''/'']';
                ftpStruc.destFormat   = '[int2str(Year) ''/'' num2str(dayNum,''%03d'') ''/'']';
                ftpStruc.fileFormat   =  {'[''COD0MGXFIN_'' int2str(Year) num2str(dayNum , ''%03i'') ''0000_01D_01D_OSB.BIA.gz'']' };
                ftpStruc.unzipFlag    = 1;
                
                filenamei = [settings.dcbDir eval(ftpStruc.destFormat) eval(ftpStruc.fileFormat{1})];
                filenamei = filenamei(1:end-3);
                
                dcbDatai = utility.readfiles.parseDcbBsxFile(filenamei,0);
            end
            
            filesRead = {filenamei};
            
        case 5
            if isempty(statCode)
                statCode = 'STFU'; % default to stanford
            end
            % stanford .mat file
            % save this one
            DfileNameFormat = ['%4d/%03d/' statCode 'MGXRAP_%4d%03d0000_01D_01D_DCB.mat'];
            
            % File location
            DpathName = settings.dcbMgexDir;
            
            % Filename
            DFileName = sprintf(DfileNameFormat, Year, dayNum,Year, dayNum);
            
            filenamei = [DpathName DFileName];
            if exist(filenamei,'file')
                
                temp = load(filenamei);
                dcbDatai = temp.dcbData;
            else
               dcbDatai = []; 
            end
            
            filesRead = {filenamei};
            
        case 6 % CODE rolling updates 
            DpathName      = settings.dcbDir;
            
            % the file name is always the same
            DFileName = 'code.bia';
    
            filenamei = [DpathName DFileName];
            
            dcbDatai = utility.readfiles.parseDcbBsxFile(filenamei,0); 
    end
    
    if ~isempty(dcbDatai) && (ddx == 1 || isempty(dcbData))
        dcbData = dcbDatai;
        fields = fieldnames(dcbData);
    elseif ~isempty(dcbDatai) 
        for idx = 1:length(fields)
            dcbData.(fields{idx}) = [dcbData.(fields{idx}); dcbDatai.(fields{idx});];
        end
    end
    
end

end