function initDcb(obj,year,doy,varargin)
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
DOWNLOAD        = res.DOWNLOAD;        % indicator to check for downloads and download

settings = obj.settings;

dcbType = 3; % 1 = CODE, 0/2 = SU, 3 = DLR
dcbType = 1; % 1 = CODE, 0/2 = SU, 3 = DLR

typeMap = [0 1 3;
           5 4 2];
   
type2 = typeMap(2,typeMap(1,:) == dcbType);

if DOWNLOAD
    [~,filenameDcb] =  navsu.readfiles.loadDcb(year,doy,settings,1,type2);

    if isempty(filenameDcb) || ~exist(filenameDcb{1},'file')
        navsu.ftp.download(6,year,doy,settings);
    end   
end

% Still need to pull SU estimates
if dcbType == 0
    [~,filenameDcb] =  loadDcb(year,doy,settings,1,5);
    
    if ~exist(filenameDcb{1},'file')
        % generate a dcb file
        
        consts = settings.multiConst;
        
        dcbData = genDcbEstTECMap(doy,year,'STFU',consts,settings);
    else
        [dcbData,filenameDcb] =  navsu.readfiles.loadDcb(year,doy,settings,0,5);
    end
    dcbType = 2;
elseif dcbType == 1
    % CODE
   dcbData = navsu.readfiles.loadDcb(year,doy,settings,0,4);
elseif dcbType == 3
    % DLR
    dcbData = navsu.readfiles.loadDcb(year,doy,settings,0,2);
end


obj.dcb = dcbData;
obj.dcb.type = dcbType;


end