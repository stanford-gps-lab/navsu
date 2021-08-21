function [ionoData,dcbData,IFileName,IFileNameFull] = loadIonex(yearList,dayList,settings,FLAG_NO_LOAD,center)
% loadIonex
% DESCRIPTION:
%   Find and parse IGS differential code bias corrections and ionospheric
%   data from .  The files to be parsed should already exist locally.
%
% INPUTS:
%  yearList            - N-length vector of years of desired outputs
%  dayList             - N-length vector of days of years of desired
%                        outputs
%  settings            - settings structure
%   .dcbDir            - Directory containing precise products- should be
%                        setup in initSettings with a config file
%
% OPTIONAL INPUTS:
%  FLAG_NO_LOAD        - True = do not parse the file, just output the name
%                        and locaiton of the local file
%  center              - 3 letter IGS AC code string
%
% OUTPUTS:
%  ionoData            - Structure containing ionospheric map from IONEX
%                        file
%  dcbData             - Structure containing parsed differential code bias
%                        information
%  IFileName           - Name of files parsed
%  IFileNameFull       - Full names of files parsed including paths
%
% See also: navsu.ftp.download, navsu.svOrbitClock, navsu.readfiles.loadDcb

if nargin < 4
    FLAG_NO_LOAD = 0;
end

if nargin < 5
    center = 'cod';
end

if length(dayList) > 1
    IFileName = {}; IFileNameFull = {};
    ionoDataTemp = []; dcbData = [];
    
    for idx = 1:length(dayList)
        [ionoDatai,~,IFileNamei,IFileNameFulli] = navsu.readfiles.loadIonex(yearList(idx),dayList(idx),settings,FLAG_NO_LOAD,center);
        
        IFileNameFull = [IFileNameFull; IFileNameFulli];
        IFileName = [IFileName; IFileNamei];
        ionoDataTemp{idx} = ionoDatai;
    end
    
    % combine the data that was just parsed
    ionoData = [];
    for idx = length(dayList):-1:1
        if ~isempty(ionoDataTemp{idx})
            if isempty(ionoData)
                % Initialize
                ionoData = ionoDataTemp{idx};
            else
                % Add the new data to the previously existing structure
                ionoDatai = ionoDataTemp{idx};
                ionoData.tecMap = cat(4,ionoData.tecMap,ionoDatai.tecMap);
                ionoData.epochs = cat(1,ionoData.epochs,ionoDatai.epochs);
                ionoData.datevecs = cat(1,ionoData.datevecs,ionoDatai.datevecs);
            end
        end
    end
    if ~isempty(ionoData)
        % Sort everything just in case
        [~,ixb] = sort(ionoData.epochs);
        ionoData.tecMap = ionoData.tecMap(:,:,:,ixb);
        ionoData.epochs = ionoData.epochs(ixb);
        ionoData.datevecs = ionoData.datevecs(ixb,:);
    end
    
else
    ionoData = [];
    dcbData = [];
        
    target_dir = fullfile(settings.dcbDir, ...
                          int2str(yearList), ...
                          num2str(dayList,'%03d'));
    
    filenamei = [center 'g' num2str(dayList,'%03d') '0.' num2str(mod(yearList,100),'%02d') 'i'];


    filenameFulli = fullfile(target_dir, filenamei);
    
    IFileName = {filenamei};
    IFileNameFull = {filenameFulli};
    
    if ~FLAG_NO_LOAD
        [dcbData,ionoData] = navsu.readfiles.readIonex(filenameFulli);
    end
end

end