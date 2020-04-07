function outDesc = convertConstIndName(inDesc,letterFlag)
% convertConstIndName
% DESCRIPTION:
%   Given the constellation index (G = 1, R = 2, E = 3, C = 4, S = 5), output 
%   the name of the constellation.  Or vice versa!
% INPUT:
%   inDesc      - Either the constellation index or the name of the
%                 constellation.
%
% OPTIONAL INPUTS:
%   letterFlag  - If true, the names of the constellations used are actually
%                 just the single letter versions, i.e. 'GPS' -> 'G'
% OUTPUT:
%   outDesc     - Either the name of the constellation or the constellation
%                 index, depending on the input. 
%
% See also: navsu.svprn.constSvnData,


consts = {'GPS','GLO','GAL','BDS','SBAS','QZSS'};

if nargin < 2
    letterFlag = false;
end

if letterFlag
    consts = {'G' 'R' 'E' 'C' 'S' 'J'};
end

if ischar(inDesc)
    outDesc = find(contains(consts,inDesc));
else 
    outDesc = consts{inDesc};
end

end