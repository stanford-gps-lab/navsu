function biasOut = findDcbElement(PRNs,consts,obs1,obs2,epochs,dcbData,statCode,noNan)
% findDcbElement
% DESCRIPTION:
%   Given DCB data parsed and placed into a structure, pull the correct
%   bias for a specific satellite, time, and observation or pair of observations 
% INPUT:
%   PRNs    - vector list of PRNs to find biases for
%   consts  - vector list of corresponding constellations- this is
%             constellation index (GPS = G = 1, GRECS = 12345)
%   obs1    - cell array of strings of RINEX 3 observation codes 
%   obs2    - cell array of strings of RINEX 3 observation codes (or 'ABS')
%   epochs  - reference times for each of the desired outputs- GPS epoch
%   dcbData - parsed DCB stucture from navsu.readfiles.loadDcb
%   statCode- optional input for receiver station code i.e. 'STFU'
%   
% OUTPUT:
%   biasOut - vector list of DCB biases in seconds
%
% See also: navsu.readfiles.loadDcb

if nargin < 7
    statCode = [];
end

if nargin < 8
    noNan = false;
end

if noNan
    biasOut = zeros(size(PRNs));
else
    biasOut = nan(size(PRNs));
end

for idx = 1:length(PRNs)
    prni = PRNs(idx);
    consti = consts(idx);
    
    obs1i = obs1{idx};
    obs2i = obs2{idx};
    
    epochi = epochs(idx);
    
    if isnan(epochi)
        epochi = inf;
    end
    epoch1 = epochi;
    epoch2 = epochi;
    
    inds = dcbData.PRNs == prni ...
         & dcbData.constInd == consti ...
         & navsu.internal.strFindCell(dcbData.obs1,obs1i,1,1) ...
         & navsu.internal.strFindCell(dcbData.obs2,obs2i,1,1) ...
         & epoch1 >= dcbData.startEpoch ...
         & epoch2 <= dcbData.endEpoch;

    if ~isempty(statCode)
        inds = inds ...
             & navsu.internal.strFindCell(dcbData.sites,statCode,1,1);
    end
    
    if any(inds)
        biasOut(idx) = dcbData.bias(find(inds, 1))*1e-9;
    end
    
end


end