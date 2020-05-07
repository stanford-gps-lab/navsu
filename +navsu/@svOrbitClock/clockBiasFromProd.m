function cbias = clockBiasFromProd(obj,prns,constInds,epochs)

% check whether to use precise clock

orbClockData = obj;

if strcmp(orbClockData.clkMode,'PRECISE')
    cbias = obj.clockInterp(prns,constInds,epochs,orbClockData.PClock);
else
    cbias = nan(size(prns));
   % put broadcast clock stuff here 
      % do broadcast stuff
    constUn = unique(constInds);
    
    % convert all epochs to week number and TOW
    [weeks,tows] = navsu.time.epochs2gps(epochs);
    
    if ismember(1,constUn)
        % gps    
        indsi = find(constInds == 1);
        pos = navsu.geo.propNavMsg(orbClockData.BEph.gps,prns(indsi),weeks(indsi),tows(indsi),'GPS');
        cbias(indsi) = pos.clock_bias-pos.TGD*0;
    end
    
    if ismember(2,constUn)
        % glo    
        indsi = find(constInds == 2);
        pos = navsu.geo.propNavMsgGlo(orbClockData.BEph.glo,prns(indsi),weeks(indsi),tows(indsi));
        cbias(indsi) = pos.clock_bias;
    end
    
    if ismember(3,constUn)
        % galileo    
        indsi = find(constInds == 3);
        pos = navsu.geo.propNavMsg(orbClockData.BEph.gal,prns(indsi),weeks(indsi),tows(indsi),'GAL');
        cbias(indsi) = pos.clock_bias-pos.TGD*0;
    end
    
     if ismember(4,constUn)
        % beidou    
        indsi = find(constInds == 4);
        pos = navsu.geo.propNavMsg(orbClockData.BEph.bds,prns(indsi),weeks(indsi),tows(indsi),'BDS');
        cbias(indsi) = pos.clock_bias-pos.TGD*0;
    end
    
end



end