
function outCodes = convertRinex3ObsCodes(inCodes,const)
% RINEX 2->3 MAPPING
% GPS
mapGPS  =  {'C1' 'C1C';
'P1' 'C1P';
'L1' 'L1C';
'D1' 'D1C';
'S1' 'S1C';
'C2' 'C2X';
'P2' 'C2W';
'L2' 'L2W';
'D2' 'D2W';
'S2' 'S2W';
'C5' 'C5X';
'L5' 'L5X';
'D5' 'D5X';
'S5' 'S5X'};

mapGLO = {'C1' 'C1C';
'P1' 'C1P';
'L1' 'L1C';
'D1' 'D1C';
'S1' 'S1C';
'C2' 'C2C';
'P2' 'C2P';
'L2' 'L2C';
'D2' 'D2C';
'S2' 'S2C'};

mapGAL = {'C1' 'C1X';
'L1' 'L1X';
'D1' 'D1X';
'S1' 'S1X';
'C5' 'C5X';
'L5' 'L5X';
'D5' 'D5X';
'S5' 'S5X';
'C7' 'C7X';
'L7' 'L7X';
'D7' 'D7X';
'S7' 'S7X';
'C8' 'C8X';
'L8' 'L8X';
'D8' 'D8X';
'S8' 'S8X';
'C6' 'C6X';
'L6' 'L6X';
'D6' 'D6X';
'S6' 'S6X'};

mapSBAS = {'C1' 'C1C';
'L1' 'L1C';
'D1' 'D1C';
'S1' 'S1C';
'C5' 'C5I';
'L5' 'L5I';
'D5' 'D5I';
'S5' 'S5I'};


switch const
    case 'G'
        mapi = mapGPS;
    case 'R'
        mapi = mapGLO;
    case 'E'
        mapi = mapGAL;
    case 'C'
        mapi = [];
    case 'S'
        mapi = mapSBAS;
    otherwise
        mapi = [];
end

outCodes = repmat({''},length(inCodes),1);
if isempty(mapi)
    return;
end

for idx = 1:size(inCodes)
    indi = find(~cellfun(@isempty, strfind(mapi(:,1),inCodes{idx})));
    
    if ~isempty(indi)
        outCodes{idx} = mapi{indi,2};
    end    
end



end

