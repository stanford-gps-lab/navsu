function initClockData(obj,years,doys,varargin)
% Download (if necessary) and load GNSS high rate (5 or 30 sec) clock data
% DESCRIPTION:
%   Initialize the precise clock structure from navsu.svOrbitClock.  This
%   will load from local files or download from IGS FTP sites if necessary.
%   
% INPUT:
%   years   - vector of years associated with the days of years when data
%             should be downloaded
%   doys    - vector of days of year for which data should be downloaded
%  
% OPTIONAL INPUTS:
%   DOWNLOAD - boolean indcating whether or not to try to download the
%             required products if they aren't found locally. Default is
%             TRUE. 
%             
% OUTPUT:
%   The object will have a brand new clock structure!  
%
% See also: navsu.svOrbitClock.clock, navsu.svOrbitClock.initOrbitData


p = inputParser;
p.addParameter('DOWNLOAD',true);

% parse the results
parse(p, varargin{:});
res = p.Results;
DOWNLOAD      = res.DOWNLOAD;         % flag to automatically download products if needed

settings = obj.settings;

% if there is no precise ephemeris already there, just make a new one
if isempty(obj.PClock)
    
    if DOWNLOAD
        [~,filenames,filenameFull] = navsu.readfiles.loadCFst(years,doys,settings,true);
        
        % build a list of IGS AC codes to download
        fileAvailable = cellfun(@exist,filenameFull);
        igsCodes = cellfun(@(x) x(1:3),filenames,'UniformOutput',false);
        
        codesDownload = unique(igsCodes(~fileAvailable));
        
        for idx =1 :length(codesDownload)
           navsu.ftp.download(2,years,doys,settings,codesDownload{idx});
        end
        
    end
    obj.PClock = navsu.readfiles.loadCFst(years,doys,settings);
else
    % This will be added to the rest of it
    
    % ADD STUFF HERE PLEASE
end




end