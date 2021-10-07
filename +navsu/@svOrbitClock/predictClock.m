function cbias = predictClock(obj,prns,constInds,epochs,latency)
 

% Initialize

cbias = nan(size(prns,1),1);

% Settings that should be pulled out
tFit = 3600*24;
orderFitClock = 1;

nProp = size(prns,1);

for idx = 1:nProp
    
    prni = prns(idx);
    consti = constInds(idx);
    epochi = epochs(idx);
    
    % Find PRN index
    indPrn = find(obj.PClock.PRNs == prni & obj.PClock.constInds == consti);
    
    % Propagation time
    epochFitEnd = min([epochi-latency max(obj.PClock.Cepochs)]);
    epochFitStart = epochFitEnd-tFit;
    
   
    
    indsFit = find(obj.PClock.Cepochs >= epochFitStart & obj.PClock.Cepochs <= epochFitEnd);
    
    epochsi = obj.PClock.Cepochs(indsFit);
    clockBase = obj.PClock.Cclk(indPrn,indsFit)';
    
    cbias(idx) = navsu.geo.polyinterp(epochsi,clockBase,orderFitClock,epochi);
end

 
 end