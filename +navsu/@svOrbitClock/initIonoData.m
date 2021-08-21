function initIonoData(obj,year,doy,varargin)
% Download (if necessary) and load ionospheric data
% DESCRIPTION:
%   Initialize the ionospheric structure from svOrbitClock
%   
% INPUT:
%   years   - vector of years associated with the days of years when data
%             should be downloaded
%   doys    - vector of days of year for which data should be downloaded
%             
% OUTPUT:
%   The object will have a brand new iono structure!  
%
% See also: navsu.svOrbitClock1

settings = obj.settings;

[~,~,~,IFileNameFull] = navsu.readfiles.loadIonex(year,doy,settings,1);

% check if we need to download anything
dlFlag = ~cellfun(@isfile, IFileNameFull);

if any(dlFlag)
    % Download necessary iono products
    navsu.ftp.download(15, year(dlFlag), doy(dlFlag), settings);
end

% now parse the files
[ionoData,~,~,filenameIono] = navsu.readfiles.loadIonex(year,doy,settings,0);

obj.iono = ionoData;

end