function initDcb(obj,year,doy,varargin)


p = inputParser;

p.addParameter('DOWNLOAD',true);   

% parse the results
parse(p, varargin{:});
res = p.Results;
DOWNLOAD        = res.DOWNLOAD;        % indicator to check for downloads and download

settings = obj.settings;

dcbType = 3; % 1 = CODE, 0/2 = SU, 3 = DLR

typeMap = [0 1 3;
           5 4 2];
   
type2 = typeMap(2,typeMap(1,:) == dcbType);

if DOWNLOAD
    [~,filenameDcb] =  loadDcb(year,doy,settings,1,type2);

    if isempty(filenameDcb) || ~exist(filenameDcb{1},'file')
        ftpHelper(6,year,doy,settings);
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
        [dcbData,filenameDcb] =  loadDcb(year,doy,settings,0,5);
    end
    dcbType = 2;
elseif dcbType == 1
    % CODE
   dcbData = loadDcb(year,doy,settings,0,4);
elseif dcbType == 3
    % DLR
    dcbData = loadDcb(year,doy,settings,0,2);
end


obj.dcb = dcbData;
obj.dcb.type = dcbType;




end