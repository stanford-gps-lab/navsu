function cbias = clock(obj,prns,constInds,epochs,varargin)
% Primary interpolation method for precise orbits
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
%   latency   - artificial latency for orbit prediction
%
% OUTPUT:
%   cbias     - N-length vector of satellite clock biases [s]
%   
% See also: navsu.svOrbitClock.initClockData, navsu.svOrbitClock.propagate
%           navsu.svOrbitClock.initClockData,

% this is mostly a wrapper for clockBiasFromProd
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