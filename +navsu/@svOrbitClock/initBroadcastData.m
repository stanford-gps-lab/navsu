function initBroadcastData(obj,years,doys,varargin)


p = inputParser;
p.addParameter('DOWNLOAD',true);   

% parse the results
parse(p, varargin{:});
res = p.Results;
DOWNLOAD        = res.DOWNLOAD;        % indicator to check for downloads and download

settings = obj.settings;

if DOWNLOAD
    [~,filenames,filenameFull] = navsu.readfiles.loadBrdc(years,doys,settings,'FLAG_NO_LOAD',true);
    
    % build a list of IGS AC codes to download
    fileAvailable = cellfun(@exist,filenameFull);
   
    navsu.ftp.download(16,years,doys,settings);
end

% Load the data!
[eph,filenames,filenameFull] = navsu.readfiles.loadBrdc(years,doys,settings);

% Put it in the object
obj.BEph = eph;

end