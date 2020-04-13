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
dlFlag = 0;
for idx = 1:length(IFileNameFull)
    if ~exist(IFileNameFull{idx},'file')
        dlFlag = 1;
    end
end
if dlFlag
    % Download necessary iono products
    navsu.ftp.download(15,year,doy,settings);
end
[ionoData,~,~,filenameIono] = navsu.readfiles.loadIonex(year,doy,settings,0);

obj.iono = ionoData;

end