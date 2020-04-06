function cbias = clock(obj,prns,constInds,epochs,varargin)


% this is mostly a wrapper for svPosFromProd
p = inputParser;
p.addParameter('latency',0);


% parse the results
parse(p, varargin{:});
res = p.Results;
latency = res.latency;
     

if ~strcmp(obj.clkMode,'PREDICT')
    cbias =  obj.clockBiasFromProd(prns,constInds,epochs);
else
    cbias = obj.predictClock(prns,constInds,epochs,latency); 
end


end