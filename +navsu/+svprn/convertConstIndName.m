function outDesc = convertConstIndName(inDesc,letterFlag)
% convertConstIndName
% DESCRIPTION:
%   Given the constellation index (G = 1, R = 2, E = 3, C = 4, J = 5, S = 6),
%   output the name of the constellation.  Or vice versa!
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

if nargin < 2
    letterFlag = false;
end

if letterFlag
    consts = {'G' 'R' 'E' 'C' 'J' 'S'};
else
    consts = {'GPS','GLO','GAL','BDS','QZSS','SBAS'};
end

if ischar(inDesc)
    outDesc = find(contains(consts, inDesc));
else
    if ~islogical(inDesc) && all(ismember(inDesc, [0 1]))
        % double array of 1's and 0's that is mimicing a logical array
        inDesc = logical(inDesc);
    end
    outDesc = consts{inDesc};
end

end