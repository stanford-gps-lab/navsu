function cbias = clock(obj,prns,constInds,epochs,varargin)


% this is mostly a wrapper for svPosFromProd
p = inputParser;

% p.addParameter('sunPos',[]);
% p.addParameter('atxData',[]);
% p.addParameter('FLAG_APC_OFFSET',true);
% p.addParameter('pPosInds',[]);
% p.addParameter('pPosPoly',[]);
% p.addParameter('dttx',0);

% parse the results
parse(p, varargin{:});
res = p.Results;
% sunPos          = res.sunPos;          % ECEF position of the sun
% atxData         = res.atxData;         % IGS ANTEX data
% FLAG_APC_OFFSET = res.FLAG_APC_OFFSET; % whether or not to add the antenna phase center offset
% pPosInds        = res.pPosInds;        
% pPosPoly        = res.pPosPoly;         
% dttx            = res.dttx;            


cbias =  clockBiasFromProd(prns,constInds,epochs,obj);






end