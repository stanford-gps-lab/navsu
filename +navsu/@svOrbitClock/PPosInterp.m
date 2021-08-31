function [posP, velP, pPosInds, pPosPoly,ROut] = PPosInterp(obj, PRNs, epochs,...
    settings, Pvel, pPosInds, pPosPoly, constInds, FLAG_APC_OFFSET,...
    atxData, sunPos, dttx)

%% PPosInterp
% Given IGS precise orbit from IGS products, this function performs either
% lagrange or polynomial interpolation of the orbit to the desired epochs
%
% Required Inputs:
%  PRNs                - N-length vector of PRN per desired interpolation
%  epochs              - N-length GPS epoch per desired interpolation
%  settings            - settings structure that contains...
%   .nPolyFit          - number of points to use for the polynomial fit
%   .pfit              - order of fit for polynomial interpolation
%   .orbitInterpMet    - 'poly' or 'lagrange' 
%
% Optional Inputs:
%  Pvel
%  pPosInds
%  pPosPoly
%  constInds           - N-length vector of constellation index per desired
%                        interpolation.  Defaults to all GPS
%  FLAG_APC_OFFSET     - flag indicating whether or not to offset the
%                        output position by the antenna phase center using 
%                        a nominal attitude model
%  atxData             - stucture containing data from IGS ATX file- can be
%                        produced by navsu.readfiles.readAtx.m
%  sunPos              - position of the sun in ECEF
%  dttx                - additional fine time offset
%
% Outputs:
%  posP               - interpolated satellite precise position
%  velP               - interpolated satellite precise velocity
%  pPosInds           - 
%  pPosPoly           -
%  ROut               - rotation matrix from ECEF to SV body frame (or vice
%                       versa)


%%

% pull some ephemeris parameters from object
Ppos = obj.PEph.position;
Pprns = obj.PEph.PRN;
Pepochs = obj.PEph.epochs;
PconstInds = obj.PEph.constellation;

% Set up interpolation constants
n     = settings.polyfit.nPolyFit; % number of nominal points to use for polynominal fit
pfit  = settings.polyfit.pfit;     % order of fit for position interpolation

% flags for velocity computation
velCalc = nargout > 1;
noVelData = nargin < 5 || isempty(Pvel) || all(isnan(Pvel), 'all');

if nargin < 6 || isempty(pPosInds)
    pPosInds = zeros(length(epochs), 2);
end

if nargin < 7 || isempty(pPosPoly)
    pPosPoly = NaN(length(epochs), 6, pfit+1);
end

if nargin < 8
    multiConst = 0;
else
    multiConst = 1;
end

if nargin < 9
    FLAG_APC_OFFSET = 0;
end

if nargin < 11 || isempty(sunPos)
    computeSunPosFlag = 1;
else
    computeSunPosFlag = 0;
end

if nargin < 12 || isempty(dttx)
    dttx = zeros(length(epochs),1);
end

% dt to use when there is no velocity data available
dt = 0.1;


% Initialize
posP = nan(size(PRNs,1),3);
velP = nan(size(PRNs,1),3);
ROut = nan(size(PRNs,1),3,3);

for idx = 1:length(PRNs)
    prn = PRNs(idx);
    epochi = epochs(idx);
    dti = dttx(idx);
    
    if prn == 0 || isnan(epochi)
        continue;
    end
    
    % retrieve orbit positions across which to interpolate
    if multiConst
        indi = Pprns == prn & PconstInds == constInds(idx);
    else
        indi = Pprns == prn;
    end
    
    if ndims(Ppos) == 3
        Pposi = squeeze(Ppos(:, :, indi));
        Pepochsi = Pepochs;
    else
        Pposi = Ppos(indi, :);
        Pepochsi = Pepochs(indi);
    end
    
    
    % get indices of interpolation start/end
    switch settings.orbitInterpMethod
        case 'poly'
            ind1 = find(Pepochsi <= epochi, 1, 'last') - n/2 + 1;
            ind2 = find(Pepochsi > epochi, 1, 'first') + n/2 - 1;
        case 'lagrange'
            ind1 = find(Pepochsi <= epochi, 1, 'last') - floor(n/2);
            ind2 = find(Pepochsi > epochi, 1, 'first') + floor(n/2) - 1;
    end
    if (isempty(ind1) && isempty(ind2)) || any(isnan(Pposi(ind1:ind2, :)), 'all') 
        continue
    end
    
    %% now interpolate orbit position and velocity
    switch settings.orbitInterpMethod
        case 'poly'
            [posP(idx, :), ~, ~, pPosPoly(prn, 1:3, :)] = navsu.geo.polyinterp(Pepochsi(ind1:ind2) - epochi, Pposi(ind1:ind2,:), pfit, dti);
        case 'lagrange'
            posP(idx, :) = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)'-epochi, Pposi(ind1:ind2,:)', dti)';
    end

    if velCalc
        if noVelData
            % No precise velocity data available... interpolate
            % position and just use velocity based on that :(
            switch settings.orbitInterpMethod
                case 'poly'
                    posPP = navsu.geo.polyinterp(Pepochsi(ind1:ind2)-epochi, Pposi(ind1:ind2,:), pfit, dt+dti);
                case 'lagrange'
                    posPP = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)'-epochi, Pposi(ind1:ind2,:)', dt+dti)';
            end

            velP(idx, :) = (posPP - posP(idx, :)) / dt;
        else
            Pveli = Pvel(indi, :);

            [velP(idx, :),~,~,pPosPoly(prn,4:6,:)] = navsu.geo.polyinterp(Pepochsi(ind1:ind2)-epochi, Pveli(ind1:ind2,:), pfit, epochi);

        end
    end

    pPosInds(prn,1) = ind1;
    pPosInds(prn,2) = ind2;
end


if FLAG_APC_OFFSET  && ~isempty(atxData)
    % Compute sun position
    typeVec = [atxData.type];
    prnVec  = [atxData.prn];
    epochStartVec = [atxData.epochStart];
    epochEndVec   = [atxData.epochEnd];
    
    offsetECEF = zeros(3,length(PRNs));
    for pdx = 1:length(PRNs)
        if computeSunPosFlag
            if pdx == 1
                lastEpoch = 0;
            end
            if lastEpoch ~= round(epochs(pdx))
                jd = navsu.time.epochs2jd(epochs(pdx));
%                 et       = cspice_str2et(['jd ' num2str(jd)]);
%                 sunposi  = cspice_spkezr( 'sun',et , 'itrf93', 'none', 'earth');
%                 sunpos = sunposi(1:3)*1000;
                
                sunpos = navsu.geo.sunVecEcef(jd)';
                lastEpoch = round(epochs(pdx));
            end
            
        else
            sunpos = sunPos(:,pdx);
        end
        
        adx = find( typeVec == constInds(pdx) & prnVec == PRNs(pdx) & epochStartVec <= epochi & epochEndVec >= epochi);
        
        if ~isempty(adx) && ~isempty(atxData(adx).apc)
            offset = (atxData(adx).apc(1,:)')*1e-3;
        else
            continue
        end
        
        svPosi = posP(pdx,:)';
        
        e = (sunpos-svPosi)./norm(sunpos-svPosi);
        k = -svPosi./norm(svPosi);
        % yhat = cross(e,k)./norm(cross(e,k));
        j = cross(k,e)/norm(cross(k,e));
        i = cross(j,k)/norm(cross(j,k));
        
        R = [i j k];
        offsetECEF(:,pdx) = R*offset;
        ROut(pdx,:,:) = R;
        
    end
    
    posP = posP + offsetECEF';
    
end


end






















