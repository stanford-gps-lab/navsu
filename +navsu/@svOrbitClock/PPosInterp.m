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
%  Ppos                - Precise position matrix from IGS products
%  Pprns               - PRNS associated with Ppos
%  Pepochs             - GPS epochs associated with Ppos
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
%  PconstInds    
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

% pull obsolete inputs from object
Ppos = obj.PEph.position;
Pprns = obj.PEph.PRN;
Pepochs = obj.PEph.epochs;
PconstInds = obj.PEph.constellation;

% Set up interpolation constants
n     = settings.polyfit.nPolyFit; % number of nominal points to use for polynominal fit
pfit  = settings.polyfit.pfit;     % order of fit for position interpolation

if nargin < 7 || isempty(Pvel)
    velCalc = 0;
    velP = [];
else
    velCalc = 1;
end

if nargin < 8 || isempty(pPosInds)
    pPosInds = zeros(length(epochs),2);
end

if nargin < 9 || isempty(pPosPoly)
    %     pPosPoly
end

if nargin < 10
    multiConst = 0;
else
    multiConst = 1;
end

if nargin < 12
    FLAG_APC_OFFSET = 0;
end

if nargin < 14 | isempty(sunPos)
    computeSunPosFlag = 1;
else
    computeSunPosFlag = 0;
end

if nargin < 15 | isempty(dttx)
    dttx = zeros(length(epochs),1);
end

% dt to use when there is no velocity data available
dt = 0.1;
if sum(sum(isnan(Pvel))) == size(Pvel,1)*size(Pvel,2)
    noVelData = 1;
else
    noVelData = 0;
end

% Initialize
posP = nan(size(PRNs,1),3);
velP = nan(size(PRNs,1),3);
ROut = nan(size(PRNs,1),3,3);
constLast = 0;
for idx = 1:length(PRNs)
    prn = PRNs(idx);
    epochi = epochs(idx);
    dti = dttx(idx);
    
    if prn == 0 || isnan(epochi)
        continue;
    end
    
    if multiConst
        if ndims(Ppos) == 3
            indi = find(Pprns == prn & PconstInds == constInds(idx));
            Pposi = squeeze(Ppos(:,:,indi));
            Pepochsi = Pepochs;
        else
            if constInds(idx) ~= constLast
                Pposc =  Ppos( PconstInds == constInds(idx),:);
                Pepochsc = Pepochs(PconstInds == constInds(idx));
                Pprnsc  = Pprns(PconstInds == constInds(idx));
            end
            constLast = constInds(idx);
            Pposi = Pposc(Pprnsc == prn,:);
            Pepochsi = Pepochsc(Pprnsc == prn);
        end
    elseif ~multiConst
        Pposi    = Ppos(Pprns == prn,:);
        Pepochsi = Pepochs(Pprns == prn);
    end
    switch settings.orbitInterpMethod
        case 'poly'
            ind1 = max(find(Pepochsi <= epochi))-n/2+1;
            ind2 = min(find(Pepochsi > epochi))+n/2-1;
        case 'lagrange'
            ind1 = max(find(Pepochsi <= epochi))-floor(n/2);
            ind2 = min(find(Pepochsi > epochi))+floor(n/2)-1;
    end
    if (isempty(ind1) && isempty(ind2)) || (sum(isnan(Pposi(ind1:ind2))) == (ind2-ind1+1)) 
        continue
    end
    
    % Check if we already have this polynomial computed
    if multiConst || pPosInds(prn,1) ~= ind1 && pPosInds(prn,2) ~= ind2
        switch settings.orbitInterpMethod
            case 'poly'
                [posP(idx,1),~,~,pPosPoly(prn,1,:)] = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,1),pfit,epochi);
                [posP(idx,2),~,~,pPosPoly(prn,2,:)] = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,2),pfit,epochi);
                [posP(idx,3),~,~,pPosPoly(prn,3,:)] = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,3),pfit,epochi);
            case 'lagrange'
                [posP(idx,1)] = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)'-epochi,Pposi(ind1:ind2,1)',epochi-epochi+dti);
                [posP(idx,2)] = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)'-epochi,Pposi(ind1:ind2,2)',epochi-epochi+dti);
                [posP(idx,3)] = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)'-epochi,Pposi(ind1:ind2,3)',epochi-epochi+dti);
        end
        
        if velCalc
            if noVelData
                % No precise velocity data available... interpolate
                % position and just use velocity based on that :(
                switch settings.orbitInterpMethod
                    case 'poly'
                        posP1 = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,1),pfit,epochi+dt);
                        posP2 = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,2),pfit,epochi+dt);
                        posP3 = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,3),pfit,epochi+dt);
                    case 'lagrange'
                        posP1 = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)',Pposi(ind1:ind2,1)',epochi+dt);
                        posP2 = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)',Pposi(ind1:ind2,2)',epochi+dt);
                        posP3 = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)',Pposi(ind1:ind2,3)',epochi+dt);
                end
                
                velP(idx,1) = (posP1-posP(idx,1))/dt;
                velP(idx,2) = (posP2-posP(idx,2))/dt;
                velP(idx,3) = (posP3-posP(idx,3))/dt;
                
            else
                Pveli =  Pvel(Pprns == prn,:);
                
                [velP(idx,1),~,~,pPosPoly(prn,4,:)] = polyinterp(Pepochsi(ind1:ind2),Pveli(ind1:ind2,1),pfit,epochi);
                [velP(idx,2),~,~,pPosPoly(prn,5,:)] = polyinterp(Pepochsi(ind1:ind2),Pveli(ind1:ind2,2),pfit,epochi);
                [velP(idx,3),~,~,pPosPoly(prn,6,:)] = polyinterp(Pepochsi(ind1:ind2),Pveli(ind1:ind2,3),pfit,epochi);
            end
            
        end
        
        pPosInds(prn,1) = ind1;
        pPosInds(prn,2) = ind2;
        
    else
        switch settings.orbitInterpMethod
            case 'poly'
                [posP(idx,1)] = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,1),pfit,epochi,0,0,squeeze(pPosPoly(prn,1,:)));
                [posP(idx,2)] = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,2),pfit,epochi,0,0,squeeze(pPosPoly(prn,2,:)));
                [posP(idx,3)] = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,3),pfit,epochi,0,0,squeeze(pPosPoly(prn,3,:)));
            case 'lagrange'
                [posP(idx,1)] = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)',Pposi(ind1:ind2,1)',epochi);
                [posP(idx,2)] = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)',Pposi(ind1:ind2,2)',epochi);
                [posP(idx,3)] = navsu.geo.lagrangeInter(Pepochsi(ind1:ind2)',Pposi(ind1:ind2,3)',epochi);
        end
        if velCalc
            if noVelData
                % No precise velocity data available... interpolate
                % position and just use velocity based on that :(
                posP1 = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,1),pfit,epochi+dt);
                posP2 = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,2),pfit,epochi+dt);
                posP3 = polyinterp(Pepochsi(ind1:ind2),Pposi(ind1:ind2,3),pfit,epochi+dt);
                
                velP(idx,1) = (posP1-posP(idx,1))/dt;
                velP(idx,2) = (posP2-posP(idx,2))/dt;
                velP(idx,3) = (posP3-posP(idx,3))/dt;
            else
                Pveli =  Pvel(Pprns == prn,:);
                
                velP(idx,1) = polyinterp(Pepochsi(ind1:ind2),Pveli(ind1:ind2,1),pfit,epochi,0,0,squeeze(pPosPoly(prn,4,:)));
                velP(idx,2) = polyinterp(Pepochsi(ind1:ind2),Pveli(ind1:ind2,2),pfit,epochi,0,0,squeeze(pPosPoly(prn,5,:)));
                velP(idx,3) = polyinterp(Pepochsi(ind1:ind2),Pveli(ind1:ind2,3),pfit,epochi,0,0,squeeze(pPosPoly(prn,6,:)));
            end
        end
    end
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






















