function [YearChange,dayChange] = ftpFileHr(YearList,dayList,ftpStruc,varargin)
% ftpFileHr
% DESCRIPTION:
%   This should just be called by navsu.ftp.download. Downloads products 
%   or data from IGS ftp sites given a day, year, information pointing to 
%   the specific product, and information about where to put the downloaded
%   products locally.  This is different from navsu.ftp.ftpFile in that
%   this function should be used specifically for IGS "high rate" products,
%   which are stored in different directories for each hour of the day.
% INPUT:
%   yearList = year corresponding to each day of desired downloads [Nx1]
%   dayList  = day of year of desired downloads [Nx1]
%   ftpStruc = MATLAB structure containing information about the structure
%              of the ftp site as well as info about the local file
%              structure. This is populated in navsu.ftp.download
%
% OUTPUT:
%   yearChange = year corresponding to each day where products were updated  
%   dayChange  = day of year when products were updated
%
% See also: navsu.ftp.download, navsu.ftp.ftpFile

ftpSite      = ['' ftpStruc.ftpSite];
sourceFormat = ftpStruc.sourceFormat;
destDir      = ftpStruc.destDir;
destFormat   = ftpStruc.destFormat;
fileFormat   = ftpStruc.fileFormat;
unzipFlag    = ftpStruc.unzipFlag;

% Years and days when files have been updated
YearChange = [];
dayChange  = [];

if nargin >= 4
    val1 = varargin{1};
end
if nargin >= 5
    val2 = varargin{2};
end
%% Download available files that do not already exist in local directory
% from ftp
mw = ftp(ftpSite);

for ddx = 1:length(dayList)
    dayNum = dayList(ddx);
    Year = YearList(ddx);
    
    jdi = navsu.time.doy2jd(Year,dayNum);
    [yri,mni,dyi] = navsu.time.jd2cal(jdi);
    [gpsWeek,gpsTow] = navsu.time.jd2gps(jdi);
    gpsDow = floor(gpsTow/86400);
    
    % Initial week of year
    [gpsWeek0,tow0] = navsu.time.jd2gps(navsu.time.cal2jd(Year,1,1));
    %     if tow0 ~= 0
    %         gpsWeek0 = gpsWeek0+1;
    %     end
    woy = gpsWeek-gpsWeek0+1;
    
    for hod = 0:23 % Hour of day
        target_dir = [destDir eval(destFormat)];
        
        ftpDir = eval(sourceFormat);
        
        % If we're not currently in the desired folder, go there.
        if ~strcmp(cd(mw),ftpDir)
            cd(mw,ftpDir);
        end
        
        if ~exist(target_dir,'dir') && strcmp(fileFormat{1},'[''*'']')
            mget(mw, '*', target_dir);
            change = 1;
        else
            %check to see what we already have
            llist = dir(target_dir);
            
            lname = {llist(:).name};
            ldate = [llist(:).datenum];
            
            change = 0;
            rlist = dir(mw);
            for i = 1:length(rlist)
                %             try
                %only download if we do not have it or the remote version is newer
                
                serverName = rlist(i).name;
                if isempty(lname)
                    have = 0;
                else
                    [have, idx] = ismember(serverName,lname);
                end
                
                % Check if this is one of the files we want
                desired = 0;
                
                for fdx = 1:length(fileFormat)
                    desiredName = upper(regexptranslate('wildcard',eval(fileFormat{fdx})));
                    
                    desired = desired || ~isempty(regexp(upper(serverName), desiredName,'once'));
                    
                    if desired
                        break
                    end
                end
                
                if (~have || rlist(i).datenum > ldate(idx)) && desired
                    mget(mw, serverName, target_dir);
                    change = 1;
                    
                    if unzipFlag
                        if strcmpi(exti,'.GZ')
                            gunzip([target_dir '\' serverName]);
                        else
                            unzip([target_dir '\' serverName]);
                        end
                    end
                end
            end
            
        end
        if change
            disp(['File(s) updated on day ' int2str(dayNum)]);
            YearChange = [YearChange; Year];
            dayChange  = [dayChange; dayNum];
        end
    end
    
end
close (mw)


end