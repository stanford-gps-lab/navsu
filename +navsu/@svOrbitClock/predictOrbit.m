function [svPos,svVel,iod,svClock] = predictOrbit(obj,prns,constInds,epochs,latency)

% Initialize
svPos = nan(size(prns,1),3);
svVel = nan(size(prns,1),3);
iod = nan(size(prns,1),1);
svClock = nan(size(prns,1),1);

% Settings that should be pulled out
tFit = 3600*6;
orderFitOrbit = 12;
% tFitClock = 3600*2;
orderFitClock = 1;

nProp = size(prns,1);

for idx = 1:nProp
    
    prni = prns(idx);
    consti = constInds(idx);
    epochi = epochs(idx);
    
    % Find PRN index
    indPrn = find(obj.PEph.PRN == prni & obj.PEph.constellation == consti);
    
    % Propagation time
    epochFitEnd = epochi-latency;
    epochFitStart = epochFitEnd-tFit;
    
    indsFit = find(obj.PEph.epochs >= epochFitStart & obj.PEph.epochs <= epochFitEnd ...
        & ~obj.PEph.Event(:,indPrn));
    
    epochsi = obj.PEph.epochs(indsFit);
    posBase = obj.PEph.position(indsFit,:,indPrn);
    clockBase = obj.PEph.clock_bias(indsFit,indPrn);
    
    % Build polynomial over the interval excluding times with events
    posi = nan(1,3);
    for ddx = 1:3
%         posi(ddx) = navsu.geo.polyinterp(epochsi,posBase(:,ddx),orderFitOrbit,epochi);
    end
    
    svPos(idx,:)= posi;
    svClock(idx) = navsu.geo.polyinterp(epochsi,clockBase,orderFitClock,epochi);
end


end