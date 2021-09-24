function offset = getAPCoffset(atxData, PRN, constInds, epochs)

%% getAPCoffset
%   
% Retrieve the antenna phase center (APC) offset for a given satellite.
% 
% Required Inputs:
%  atxData             - stucture containing data from IGS ATX file- can be
%                        produced by navsu.readfiles.readAtx.m
%  PRN                 - N-length vector of PRN
%  constInds           - N-length vector of constellation indices
%  epochs              - epochs at which to retrieve APC data. Can be
%                        scalar or N-length
% 
% Outputs:
%  offset              - 3 x N matrix of APC offset in satellite coordinate
%                        system

if numel(epochs) == 1
    epochs  = epochs * size(PRN);
end
if numel(constInds) == 1
    constInds = constInds * ones(size(PRN));
end

% retrieve atx data, collect in vectors
typeVec = [atxData.type];
svnVec  = [atxData.svn];
epochStartVec = [atxData.epochStart];
epochEndVec   = [atxData.epochEnd];

% preallocate result vec
offset = NaN(3, length(PRN));

for pdx = 1:length(PRN)
    
    
    epochi = epochs(pdx);
    
    % convert PRN to SV number
    svni = navsu.svprn.prn2svn(PRN(pdx), navsu.time.epochs2jd(epochi), constInds(pdx));
    
    adx = find( typeVec == constInds(pdx) ...
              & svnVec == svni ...
              & epochStartVec <= epochi ...
              & epochEndVec >= epochi );
    
    if ~isempty(adx) && ~isempty(atxData(adx).apc)
        offset(:, pdx) = (atxData(adx).apc(1,:))*1e-3;
    else
        continue
    end
    
    
end

end

