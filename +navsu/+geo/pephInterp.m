function [posP, velP, pPosInds, pPosPoly, ROut] = pephInterp(peph, prns, ...
    constInds, epochs, varargin)
%% pephInterp
% Given precise orbit from IGS products, this function performs either
% lagrange or polynomial interpolation of the orbit to the desired epochs
%   
% [posP, velP, pPosInds, pPosPoly, ROut] = pephInterp(peph, prns, ...
%                                              constInds, epochs, varargin)
%
% Required Inputs:
%  peph               - precise ephemeris struct. See navsu.readfiles.loadPEph
%  prns               - N-length vector of PRN per desired interpolation
%  constInds          - N-length vector of constellation index per desired
%                       interpolation.
%  epochs             - N-length GPS epoch per desired interpolation
%
% Optional Inputs (these need to be passed as name-value pairs):
%  FLAG_APC_OFFSET    - flag indicating whether or not to offset the
%                       output position by the antenna phase center using
%                       a nominal attitude model. Default is "false".
%  atxData            - stucture containing data from IGS ATX file- can be
%                       produced by navsu.readfiles.readAtx.m
%  sunPos             - position of the sun in ECEF. Used with atxData.
%  orbitInterpMethod  - orbit interpolation method. 'lagrange' (default) or
%                       'poly'
%  nPolyFit           - Number of points to use for the orbit
%                       interpolation. Default is 12.
%  pfit               - Order of polynomial to be used in 'poly'
%                       interpolation. Default is 8.
%  dttx               - additional fine time offset. Defaults to 0.
%  pPosInds           - Indices of data point used for interpolation.
%  pPosPoly           - Polynomial interpolation points from 'poly'.
%
% Outputs:
%  posP               - interpolated satellite precise position
%  velP               - interpolated satellite precise velocity
%  pPosInds           - Indices of data point used for interpolation.
%  pPosPoly           - Polynomial interpolation points from 'poly'.
%  ROut               - rotation matrix from ECEF to SV body frame (or vice
%                       versa)

%% Parse input parameters
p = inputParser;
p.addParameter('orbitInterpMethod', 'lagrange'); % 'lagrange' or 'poly'
p.addParameter('nPolyFit', 12);   % Number of points to use for the fit
p.addParameter('pfit', 8);        % Order of fit for polynomial interpolation
p.addParameter('FLAG_APC_OFFSET', false);  % Whether to offset the satellite
% Position to the antenna phase center
p.addParameter('sunPos', []);  % ECEF sun position
p.addParameter('dttx', zeros(size(prns))); % Additional fine time offset
p.addParameter('atxData', []);
p.addParameter('pPosInds', []);
p.addParameter('pPosPoly', []);

% parse the results
parse(p, varargin{:});
res = p.Results;
orbitInterpMethod = res.orbitInterpMethod;
n = res.nPolyFit;
pfit     = res.pfit;
FLAG_APC_OFFSET = res.FLAG_APC_OFFSET;
sunPos = res.sunPos;
dttx = res.dttx;
pPosInds = res.pPosInds;
pPosPoly = res.pPosPoly;
atxData = res.atxData;
%%

if FLAG_APC_OFFSET && isempty(atxData)
    error('Antenna phase center (.atx) data required to do CoM to APC offset')
end


% Pull things out of the precise ephemeris container
Ppos = peph.position;
Pvel = peph.velocity;
Pepochs = peph.epochs;
Pprns = peph.PRN;
PconstInds = peph.constellation;


% If any more than just the position is requested, compute the velocity
velCalc = nargout > 1;
noVelData = isempty(Pvel) || all(isnan(Pvel), 'all');

% dt to use when there is no velocity data available
dt = 0.01;

%%
% Initialize
posP = nan(size(prns,1), 3);
velP = nan(size(prns,1), 3);
ROut = nan(size(prns,1), 3, 3);
if isempty(pPosPoly)
    pPosPoly = NaN(max(prns, [], 'omitnan'), 6, pfit+1);
end
if isempty(pPosInds)
    pPosInds = NaN(max(prns, [], 'omitnan'), 2);
end


for idx = 1:length(prns)
    prn = prns(idx);
    epochi = epochs(idx);
    dti = dttx(idx);
    
    if prn == 0 || isnan(epochi)
        continue;
    end
    
    % retrieve orbit positions across which to interpolate
    indi = Pprns == prn & PconstInds == constInds(idx);
    
    if ~any(indi)
        continue
    end
    
    if ndims(Ppos) == 3
        Pposi = squeeze(Ppos(:, :, indi));
        Pepochsi = Pepochs;
    else
        Pposi = Ppos(indi, :);
        Pepochsi = Pepochs(indi);
    end
    
    
    % get indices of interpolation start/end
    switch orbitInterpMethod
        case 'poly'
            ind1 = find(Pepochsi <= epochi, 1, 'last') - n/2 + 1;
            ind2 = find(Pepochsi > epochi, 1, 'first') + n/2 - 1;
        case 'lagrange'
            ind1 = find(Pepochsi <= epochi, 1, 'last') - floor(n/2);
            ind2 = find(Pepochsi > epochi, 1, 'first') + floor(n/2) - 1;
    end
    
    % make sure ind1, ind2 are within legal limits
    ind1 = max(ind1, 1);
    ind2 = min(ind2, size(Pposi, 1));
    if isempty(ind1) || isempty(ind2) || any(isnan(Pposi(ind1:ind2, :)), 'all') 
        continue
    end
    
    %% now interpolate orbit position and velocity
    switch orbitInterpMethod
        case 'poly'
            [posP(idx, :), ~, ~, pPosPoly(prn, 1:3, :)] = ...
                navsu.geo.polyinterp(Pepochsi(ind1:ind2) - epochi, ...
                                     Pposi(ind1:ind2, :), pfit, dti);
        case 'lagrange'
            posP(idx, :) = ...
                navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)' - epochi, ...
                                        Pposi(ind1:ind2, :)', dti)';
    end

    if velCalc
        if noVelData
            % No precise velocity data available... interpolate
            % position and just use velocity based on that :(
            switch orbitInterpMethod
                case 'poly'
                    posPP = navsu.geo.polyinterp(Pepochsi(ind1:ind2) - epochi, ...
                                                 Pposi(ind1:ind2, :), pfit, dt+dti);
                case 'lagrange'
                    posPP = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)' - epochi, ...
                                                    Pposi(ind1:ind2, :)', dt+dti)';
            end

            velP(idx, :) = (posPP - posP(idx, :)) / dt;
        else
            
            if ndims(Pvel) == 3
                Pveli = squeeze(Pvel(:, :, indi));
            else
                Pveli = Pvel(indi, :);
            end

            switch orbitInterpMethod
                case 'poly'
                    [velP(idx, :),~,~,pPosPoly(prn, 4:6, :)] = ...
                        navsu.geo.polyinterp(Pepochsi(ind1:ind2) - epochi, ...
                                             Pveli(ind1:ind2, :), pfit, dti);
                case 'lagrange'
                    velP(idx, :) = ...
                        navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)' - epochi, ...
                                                Pveli(ind1:ind2, :)', dti)';
            end
        end
    end

    pPosInds(prn, 1) = ind1;
    pPosInds(prn, 2) = ind2;
end


if FLAG_APC_OFFSET  && ~isempty(atxData)
    
    % get offset vector for each satellite
    offset = navsu.ppp.getAPCoffset(atxData, prns, constInds, epochs);
    
    % Compute rotation matrix per satellite
    R = navsu.geo.svLocalFrame(posP, epochs, sunPos);
    % save output
    ROut = permute(R, [3, 1, 2]);

    % now update sat positions
    offsetECEF = zeros(3,length(prns));
    for pdx = 1:length(prns)
        offsetECEF(:, pdx) = R(:, :, pdx) * offset(:, pdx);        
    end
    % update satellite positions
    posP = posP + offsetECEF';
    
end


end