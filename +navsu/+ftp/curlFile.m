function [yearChange,dayChange] = curlFile(yearList,dayList,ftpStruc,netrcFile,cookieFile,varargin)
% curlFile
% DESCRIPTION:
%   This should just be called by navsu.ftp.download. Downloads products
%   or data from IGS ftp sites given a day, year, information pointing to
%   the specific product, and information about where to put the downloaded
%   products locally.
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
% See also: navsu.ftp.download, navsu.ftp.ftpFileHr

ftpSite      = ['' ftpStruc.ftpSite];
sourceFormat = ftpStruc.sourceFormat;
destDir      = ftpStruc.destDir;
destFormat   = ftpStruc.destFormat;
fileFormat   = ftpStruc.fileFormat;
unzipFlag    = ftpStruc.unzipFlag;

% Years and days when files have been updated
yearChange = [];
dayChange  = [];

if nargin >= 6
    val1 = varargin{1};
end
if nargin >= 7
    val2 = varargin{2};
end
%% Download available files that do not already exist in local directory
% from ftp
% mw = ftp(ftpSite);

for ddx = 1:length(dayList)
    dayNum = dayList(ddx);
    Year = yearList(ddx);
    
    jdi = navsu.time.doy2jd(Year,dayNum);
    [yri,mni,dyi] = navsu.time.jd2cal(jdi);
    [gpsWeek,gpsTow] = navsu.time.jd2gps(jdi);
    gpsDow = floor(gpsTow/86400);
    
    % Initial week of year
    [gpsWeek0,tow0] = navsu.time.jd2gps(navsu.time.cal2jd(Year,1,1));
    woy = gpsWeek-gpsWeek0+1;
    
    target_dir = [destDir eval(destFormat)];
    
    ftpDir = eval(sourceFormat);
    
    % If we're not currently in the desired folder, go there.
    sourcePath = [ftpSite ftpDir];
    
    %check to see what we already have
    llist = dir(target_dir);
    
    lname = {llist(:).name};
    %     ldate = [llist(:).datenum];
    
    if matches(fileFormat{1},'[''*'']')
        % download the entire directory
        change = 1;
        navsu.ftp.curlDownloadDirectory(sourcePath,target_dir,netrcFile,cookieFile)
        
    else
        change = 0;
        rlist = navsu.ftp.curlGetDirectoryContents(sourcePath,netrcFile,cookieFile);
        
        for i = 1:length(rlist)
            %only download if we do not have it or the remote version is newer
            
            serverName = rlist{i};
            if isempty(lname)
                have = 0;
            else
                [have, idx] = ismember(serverName,lname);
            end
            
            % Check if this is one of the files we want
            desired = 0;
            
            for fdx = 1:length(fileFormat)
                desiredName = regexptranslate('wildcard',eval(fileFormat{fdx}));
                
                desired = desired || ~isempty(regexp(serverName, desiredName,'once'));
                
                if desired
                    break
                end
            end
            
            if ~have && desired
                %                 mget(mw, serverName, target_dir);
                
                navsu.ftp.curlDownloadSingleFile([sourcePath serverName],target_dir,netrcFile,cookieFile);
                change = 1;
                
                if unzipFlag
                    navsu.readfiles.unzipFile(fullfile(target_dir,  serverName));
                end
            end
            
        end
    end
    
    %     end
    if change
        disp(['File(s) updated on day ' int2str(dayNum)]);
        yearChange = [yearChange; Year];
        dayChange  = [dayChange; dayNum];
    end
    
end
% close(mw);


end