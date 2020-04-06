function initIonoData(obj,year,doy,varargin)

settings = obj.settings;

[~,~,~,IFileNameFull] = utility.readfiles.loadIonex(year,doy,settings,1);
dlFlag = 0;
for idx = 1:length(IFileNameFull)
    if ~exist(IFileNameFull{idx},'file')
        dlFlag = 1;
    end
end
if dlFlag
    % Download necessary iono products
    utility.ftp.download(15,year,doy,settings);
end
[ionoData,~,~,filenameIono] = utility.readfiles.loadIonex(year,doy,settings,0);


obj.iono = ionoData;

end