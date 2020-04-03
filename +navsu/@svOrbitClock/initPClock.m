function initPClock(obj,years,doys,varargin)

p = inputParser;
% 
% p.addParameter('FLAG_NO_LOAD',false);
% p.addParameter('atxData',[]);
p.addParameter('DOWNLOAD',true);

% parse the results
parse(p, varargin{:});
res = p.Results;
% FLAG_NO_LOAD = res.FLAG_NO_LOAD;       % Flag to not actually load the data- should probably not be used here
% atxData      = res.atxData;            % IGS ANTEX data
DOWNLOAD      = res.DOWNLOAD;         % flag to automatically download products if needed

settings = obj.settings;

% if there is no precise ephemeris already there, just make a new one
if isempty(obj.PClock)
    
    if DOWNLOAD
        [~,filenames,filenameFull] = loadCFst(years,doys,settings,true);
        
        % build a list of IGS AC codes to download
        fileAvailable = cellfun(@exist,filenameFull);
        igsCodes = cellfun(@(x) x(1:3),filenames,'UniformOutput',false);
        
        codesDownload = unique(igsCodes(~fileAvailable));
        
        for idx =1 :length(codesDownload)
           ftpHelper(2,years,doys,settings,codesDownload{idx});
        end
        
    end
    PClock = loadCFst(years,doys,settings);
    obj.PClock = PClock;
else
    % This will be added to the rest of it
    
    % ADD STUFF HERE PLEASE
end




end