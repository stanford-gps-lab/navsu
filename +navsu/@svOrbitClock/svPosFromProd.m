function [svPos,svVel,iod,sigOrbit] = svPosFromProd(obj, prns, epochs, settings, ...
    pPosInds, pPosPoly, constInds, FLAG_APC_OFFSET, atxData, sunPos, dttx)

if nargin < 8
    FLAG_APC_OFFSET = 0;
end

if nargin < 9
    atxData = [];
end

if nargin < 10
    sunPos = [];
end

if nargin < 11
    dttx = [];
end

iod = [];

if strcmp(obj.orbMode,'PRECISE')
    % this should be adjusted to allow for all of the inputs
    [svPos, svVel] = obj.PPosInterp(prns, epochs, ...
        settings, NaN, pPosInds, pPosPoly, constInds, FLAG_APC_OFFSET, ...
        atxData, sunPos, dttx);
    
    % Just setting a standard orbit sigma here :)
    sigOrbit = 0.03*ones(size(svPos,1),1);
        
else
    
    % propNavMsg handles the various constellations
    pos = navsu.geo.propNavMsg(obj.BEph, prns, constInds, epochs);
    
    svPos       = [pos.x pos.y pos.z];
    svVel       = [pos.x_dot pos.y_dot pos.z_dot];
    iod         = pos.IODC;
    iod(constInds == 2) = NaN; % GLONASS IODC is not valid
    sigOrbit    = pos.accuracy;
    sigOrbit(constInds == 2) = 10; % GLONASS does not give accuracy
    
%      warning('sigma not set here yet :/- please use URA or something')
%      sigOrbit = 10*ones(size(svPos,1),1);
    
end


end