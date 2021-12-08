function Tiono = klobuchar(ionoCorr, epoch, lat, lon, az, el)
%Copmutes approximate Ionospheric delay of pseudoranges according to the
%Klobuchar model.
%   Tiono = klobuchar(ionoCorr, epoch, lat, lon, az, el)
%   
%   Function to compute the single frequency iono delay in seconds.
%   Based on Figure 20-4 on page 127-129 of IS_GPS-200H.
%   
%   @params:
%   ionoCorr- ionospheric correction parameters broadcasted by GPS, 1x8
%             vector
%   epoch   - GPS epoch (week*86400+secofweek) or second of week or sec of
%             day
%   lat     - receiver latitude in rad
%   lon     - receiver longitude in rad
%   az      - vector of satellite azimuths in rad
%   el      - vector of satellite elevation in rad
%   
%   @outputs:
%   Tiono   - signal delay due to ionosphere in seconds

if any(~isfinite(ionoCorr))
    % default output at max input dimensions
    Tiono = NaN(size(epoch .* lat .* lon .* az .* el)); %#ok
    return
end

% make sure I have the right number of parameters
if numel(ionoCorr) ~= 8
    error('Need 8 iono correction parameters.');
end
% make sure it's a column vector
if isrow(ionoCorr)
    ionoCorr = ionoCorr';
end

% retrieve iono parameters
alpha = ionoCorr(1:4);
beta = ionoCorr(5:8);

% convert azimuth, elevation to semi-circles
A = az/pi;
E = el/pi;

% get central angle between user position and earch projection of
% ionospheric intersection point (IIP)
psi = 0.0137 ./ (E + 0.11) - 0.022;

% get geodetic latitude and longitude of earth projection of IIPs
phi = max(min(lat/pi + psi .* cos(A), 0.416), -0.416);
lambda = lon/pi + psi .* sin(A) ./ cos(phi);

% get geomagnetic latitude of earth projection of IIPs
phi_m = phi + 0.064 * cos(lambda - 1.617);

% amplitude of cosine function
AMP = reshape(max(phi_m(:).^(1:4) * alpha, 0), size(phi_m));

% get period of cosine function
PER = reshape(max(phi_m(:).^(1:4) * beta, 72000), size(phi_m));

% get value of where to evaluate cosine function
t = mod(4.32e4 * lambda + epoch, 86400);
x = 2*pi*(t - 50400) ./ PER;

% get obliquity factor
F = 1 + 16*(0.53-E).^3;

% compute iono delay estimate
Tiono = F .* (5e-9 + AMP .* (1 - x.^2/2 + x.^4/24));
Tiono(abs(x) >= 1.57) = F(abs(x) >= 1.57) * 5e-9;

% limit to elevation > 0
Tiono(el <= 0) = NaN;

end

