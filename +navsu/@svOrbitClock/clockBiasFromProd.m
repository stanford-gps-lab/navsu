function cbias = clockBiasFromProd(obj, prns, constInds, epochs)

% check whether to use precise clock

if strcmp(obj.clkMode,'PRECISE')
    cbias = obj.clockInterp(prns, constInds, epochs);
else
    % propagate navigation broadcast
    % leverage propNavMsg here
    pos = navsu.geo.propNavMsg(obj.BEph, prns, constInds, epochs);
    cbias = pos.clock_bias;
    
end

end