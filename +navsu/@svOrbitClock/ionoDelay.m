
function [tecs, delays, tecSlant] = ionoDelay(obj, epoch, llh, varargin)

% Inputs
% epoch   | GPS epoch
% llh     | user lat-long-height in degrees and meters

% Outputs
% tecs    | total electron content for each LOS
% delays  | iono code delay in meters

p = inputParser;
p.addParameter('az', []);     % azimuth from each satellite [rad]
p.addParameter('el', []);     % elevation angle from each satellite [rad]
p.addParameter('satPos', []); % ECEF satellite position from each satellite
%                               - must have this if az/el unavailable
p.addParameter('freqs', NaN);  % frequency of each signal- required to output 
%                               actual delay rather than just TEC


% parse the results
parse(p, varargin{:});
res = p.Results;
az  = res.az;
el  = res.el;
satPos = res.satPos;
freqs  = res.freqs;

% might need to setup az and el
if isempty(az) || isempty(el) 
    if isempty(satPos)
        error('Need at least satellite position to compute TEC')
    else
        % Compute az and el 
        usrPos = llh2xyz(llh);
        [el,az] = pos2elaz(usrPos,satPos);
    end
end

epoch = epoch(:);
if length(epoch) == 1 && length(el) > 1
    epoch = epoch*ones(size(el));
end

% Compute the TEC for each LOS
tecMap   = squeeze(obj.iono.tecMap);
tecs     = zeros(size(el));
obliq    = zeros(size(el));

% build TEC for each LOS
for sdx = 1:length(el)
    [latPpi, lonPpi, obliq(sdx)] = navsu.ppp.models.ionoPiercePoint( ...
        llh(1)*pi/180, llh(2)*pi/180, az(sdx), el(sdx));
    
    if isnan(latPpi)
        continue;
    end
    latPpi = latPpi*180/pi; lonPpi = lonPpi*180/pi;
    
    indsEpochs = [find(obj.iono.epochs <= epoch(sdx), 1, 'last') ...
                  find(obj.iono.epochs > epoch(sdx), 1)];
    
    tecs(sdx) = -interpn(obj.iono.latVec,...
        obj.iono.lonVec, obj.iono.epochs(indsEpochs),...
        tecMap(:,:,indsEpochs), latPpi, lonPpi, epoch(sdx));
    
end

% if we can, swap from TEC to delay (in meters)
if nargout > 1 % only do this if requested
    
    % slant TEC
    tecSlant = obliq .* tecs;

    if ~isempty(freqs)
        % Compute the delay
        delays = tecSlant*40.3*10^15./freqs.^2;
    end
end
    

end



















