function [svPos, svVel, iod, sigma] = propagate(obj, prns, constInds, ...
    epochs, varargin)
% Primary interpolation method for orbits
% DESCRIPTION:
%   Used to interpolate (typically lagrange interpolation) precise orbits
%   that have been loaded into the svOrbitClock object using
%   navsu.svOrbitClock.initOrbitData.  
% INPUT:
%   prns      - N-length vector list of PRNs to be interpolated
%   constInds - N-length vector list of constellation indices associated
%               the prns vector (1 = GPS, 2 = GLONASS, 3 = GAL, 4 = BDS, 5
%               = SBAS)
%   epochs    - N-length vector list of GPS epochs to propagate to (GPS
%               epoch is seconds since GPS time start- see
%               navsu.time.gps2epochs)
%
% OPTIONAL INPUTS: 
%   sunPos    - ECEF sun position - can be input for speed
%   atxData   - antenna phase center structure- this is nominally just
%               whatever is already loaded into the object- see
%               navsu.svOrbitClock.initAtxData
%   FLAG_APC_OFFSET - flag indicating whether or not to offset the SV
%               from the center of mass (what is in the .sp3 file) to the
%               antenna phase center
%   pPosInds  - helper interpolation indices    (prob should be removed)
%   pPosPoly  - helper interpolation polynomial (prob should be removed)
%   dttx      - time delta                      (prob should be removed)
%   latency   - artificial latency for orbit prediction
%
% OUTPUT:
%   svPos     - [Nx3] matrix of ECEF satellite positions- can be either
%               CoM or APC depending on input [m]
%   svVel     - [Nx3] matrix of ECEF satellite velocities [m/s]  
%   iod       - issue of data (used for broadcast navigation data)
%   
% See also: navsu.svOrbitClock.initOrbitData, navsu.svOrbitClock.clock
%           navsu.svOrbitClock.initClockData,

% this is mostly a wrapper for the orbit propagation
p = inputParser;

p.addParameter('sunPos',            []);
p.addParameter('FLAG_APC_OFFSET',   true);
p.addParameter('pPosInds',          []);
p.addParameter('pPosPoly',          []);
p.addParameter('dttx',              []);
p.addParameter('latency',           0);

% parse the results
parse(p, varargin{:});
res = p.Results;
latency         = res.latency;
pPosInds        = res.pPosInds;        
pPosPoly        = res.pPosPoly;         
FLAG_APC_OFFSET = res.FLAG_APC_OFFSET; % add the antenna phase center offset?
sunPos          = res.sunPos;          % ECEF position of the sun
dttx            = res.dttx;            


if strcmp(obj.orbMode, 'PREDICT')
    [svPos, svVel, iod] = obj.predictOrbit(prns, constInds, epochs, latency);
    sigma = NaN;
    
elseif strcmp(obj.orbMode, 'PRECISE')
    % Precise orbit interpolation
    [svPos, svVel] = obj.PPosInterp(prns, constInds, epochs, ...
        pPosInds, pPosPoly, FLAG_APC_OFFSET, sunPos, dttx);
    
    iod = []; % not applicable for precise orbits
    
    % Just setting a standard orbit sigma here  - IGS is very precise
    sigma = 0.03*ones(size(svPos,1),1);
    
else
    % Use broadcast file
    pos = navsu.geo.propNavMsg(obj.BEph, prns, constInds, epochs);
    
    svPos       = [pos.x pos.y pos.z];
    svVel       = [pos.x_dot pos.y_dot pos.z_dot];
    iod         = pos.IODC;
    iod(constInds == 2) = NaN; % GLONASS IODC is not valid
    sigma       = pos.accuracy;
    sigma(constInds == 2) = 10; % GLONASS does not give accuracy

end


end