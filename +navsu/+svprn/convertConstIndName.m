function outDesc = convertConstIndName(inDesc,letterFlag)

consts = {'GPS','GLO','GAL','BDS','SBAS','QZSS'};

if nargin == 2
    consts = {'G' 'R' 'E' 'C' 'S' 'J'};
end

if ischar(inDesc)
    outDesc = find(contains(consts,inDesc));
else 
    outDesc = consts{inDesc};
end

end