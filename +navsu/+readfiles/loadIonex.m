function [ionoData,dcbData,IFileName,IFileNameFull] = loadIonex(Year,dayNum,settings,FLAG_NO_LOAD,center)
% Load IONEX data

if nargin < 4
    FLAG_NO_LOAD = 0;
end

if nargin < 5
    center = 'cod';
end

if length(dayNum) > 1
    IFileName = {}; IFileNameFull = {};
    
    for idx = 1:length(dayNum)
        [ionoDatai,IFileNamei,IFileNameFulli] = navsu.readfiles.loadIonex(Year(idx),dayNum(idx),settings,FLAG_NO_LOAD,center);
    end
    
    % do something useful with all this
    
else
    ionoData = [];
    dcbData = [];
    
    destDir      = settings.dcbDir;
    destPath   = [int2str(Year) '\' num2str(dayNum,'%03d') '\'];
   
    filenamei = [center 'g' num2str(dayNum,'%03d') '0.' num2str(mod(Year,100),'%02d') 'i'];
    target_dir = [destDir destPath];

    filenameFulli = [target_dir filenamei];
    
    IFileName = {filenamei};
    IFileNameFull = {filenameFulli};
    
    if ~FLAG_NO_LOAD
        [dcbData,ionoData] = navsu.readfiles.readIONEX(filenameFulli);
    end
end

end