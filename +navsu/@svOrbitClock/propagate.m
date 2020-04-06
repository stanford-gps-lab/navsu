function [svPos,svVel,iod,svClock] = propagate(obj,prns,constInds,epochs,varargin)


% this is mostly a wrapper for svPosFromProd
p = inputParser;

p.addParameter('sunPos',[]);
p.addParameter('atxData',obj.atx);
p.addParameter('FLAG_APC_OFFSET',true);
p.addParameter('pPosInds',[]);
p.addParameter('pPosPoly',[]);
p.addParameter('dttx',[]);
p.addParameter('latency',0);

% parse the results
parse(p, varargin{:});
res = p.Results;
sunPos          = res.sunPos;          % ECEF position of the sun
atxData         = res.atxData;         % IGS ANTEX data
FLAG_APC_OFFSET = res.FLAG_APC_OFFSET; % whether or not to add the antenna phase center offset
pPosInds        = res.pPosInds;        
pPosPoly        = res.pPosPoly;         
dttx            = res.dttx;            
latency         = res.latency;

settings = obj.settings;

if ~strcmp(obj.orbMode,'PREDICT')
    [svPos,svVel,iod] = obj.svPosFromProd(prns, epochs,settings,...
        pPosInds,pPosPoly,constInds,FLAG_APC_OFFSET,atxData,sunPos,dttx);
else
    [svPos,svVel,iod,svClock] = obj.predictOrbit(prns,constInds,epochs,latency);
end


end