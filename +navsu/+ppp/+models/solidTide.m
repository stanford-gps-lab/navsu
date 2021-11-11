function [posOffset,rangeOffset] = solidTide(epoch,pos,varargin)

% Parse optional inputs
p = inputParser;
p.addParameter('sunPos',  []); % where the sun at (ECEF, meters)
p.addParameter('moonPos', []); % where the moon at (ECEF, meters)
p.addParameter('svPos',   []); % Nx3 where the satellites at (ECEF, meters)

parse(p, varargin{:});
res = p.Results;
sunPos  = res.sunPos;
moonPos = res.moonPos;
svPos   = res.svPos;

% ensure pos is column vector
if isrow(pos)
    pos = pos';
end

%% If sun and moon position are not available, need to find them
if isempty(sunPos)
    sunPos  = navsu.geo.sunVecEcef(navsu.time.epochs2jd(epoch));
end
if isempty(moonPos)
    moonPos = navsu.geo.moonVecEcef(navsu.time.epochs2jd(epoch));
end

%% need time in year, month, day, hour of day
[yr,mn,dy,hour,min,sec] = navsu.time.epochs2cal(epoch);
hour = hour+min/60+sec/3600;

% Solid tides
posOffset = navsu.ppp.models.iers.dehanttideinel(pos,yr,mn,dy,hour,sunPos,moonPos);

%% If satellite positions are available and the output is requested, also 
% include the offsets per range
if nargout > 1
    if isempty(svPos)
        error('Satellite positions required to output solid tide range offset');
    end
    los = pos-svPos';
    los = los ./ (sqrt(sum(los.^2, 1)));
    % dot product for range offsets
    rangeOffset = sum(los.*posOffset, 1)';    
end

end




















