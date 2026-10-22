function rGeo = geocentricRadius(pos)
%rGeo = geocentricRadius(position)
%   Calculates the geocentric radius at a given position.
%   Position must be given by an N x 3 vector in ECEF coordinates in meter.
%   returns an N x 1 column vector of geocentric radi.

% make sure pos is row vector
if size(pos, 2) ~= 3
    pos = pos';
end

% Constants (some of these could be changed to inputs at a later date)
rEarth = navsu.constants.rEarthWgs84;
e = navsu.constants.eEarthWgs84;

% compute geocentric radius
llhi = navsu.geo.xyz2llh(pos);
latRad = llhi(:, 1)*pi/180;
rGeo = rEarth ./ sqrt(1 - (e * sin(latRad)).^2) ...
    .* sqrt(cos(latRad).^2 + (1 - e^2)^2 * sin(latRad).^2); % from (2.137)
end

