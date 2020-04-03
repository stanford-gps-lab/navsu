function [Obs_types, pos_M, ifound_types, interval, sysId, antoff, antmod,...
    mixedFlag, rxmod,rinexVer] = rinexParseHeader(file)

% SYNTAX:
%   [Obs_types, pos_M, ifound_types, interval, sysId, antoff, antmod] = RINEX_parse_hdr(file);
%
% INPUT:
%   file = pointer to RINEX observation file
%
% OUTPUT:
%   Obs_types = cell of strings containing observation types
%               RINEX v2.xx --> e.g. L1C1P1...
%               RINEX v3.xx --> e.g. C1CL1CD1C...
%   pos_M = master station approximate position
%   ifound_types = boolean variable to check the correct acquisition of basic information
%   interval = observation interval in seconds
%   sysId = cell-array containing one-letter identifiers for constellations
%   antoff = antenna offset [m]
%   antmod = antenna model [string]
%
% DESCRIPTION:
%   RINEX observation file header analysis.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) Kai Borre
% Kai Borre 09-23-97
%
% Adapted by Mirko Reguzzoni, Eugenio Realini, 2009
% Portions of code contributed by Damiano Triglione, 2012
%----------------------------------------------------------------------------------------------

ifound_types = 0;
Obs_types = cell(0,0);
sysId = cell(0,0);
pos_M = [];
interval = 0;
antmod = '';
rxmod = '';

%parse first line
line = fgetl(file);

mixedFlag = ~isempty(strfind(line,'MIXED'));

rinexVer = 0;

%constellation counter for RINEX v3.xx
c = 1;

%check if the end of the header or the end of the file has been reached
while isempty(strfind(line,'END OF HEADER')) && ischar(line)
    %NOTE: ischar is better than checking if line is the number -1.
    
    answer = strfind(line,'# / TYPES OF OBSERV'); %RINEX v2.xx
    if ~isempty(answer)
        rinexVer = 2;
        Obs_types{1} = [];
        nObs = sscanf(line(1:6),'%d');
        nLinObs = ceil(nObs/9);
        for i = 1 : nLinObs
            if (i > 1)
                line = fgetl(file);
            end
            n = min(nObs,9);
            for k = 1 : n
                ot = sscanf(line(k*6+1:k*6+6),'%s');
                Obs_types{1} = [Obs_types{1} ot];
            end
            nObs = nObs - 9;
        end
        
        ifound_types = 1;
    end
    
    answer = strfind(line,'SYS / # / OBS TYPES'); %RINEX v3.xx
    if ~isempty(answer)
        rinexVer = 3;
        sysId{c} = sscanf(line(1),'%s');
        nObs = sscanf(line(2:6),'%d');
        Obs_types.(sysId{c}) = [];
        nLinObs = ceil(nObs/13);
        for i = 1 : nLinObs
            if (i > 1)
                line = fgetl(file);
            end
            n = min(nObs,13);
            for k = 0 : n-1
                ot = sscanf(line(6+k*4+1:6+k*4+4),'%s');
                Obs_types.(sysId{c}) = [Obs_types.(sysId{c}) ot];
            end
            nObs = nObs - 13;
        end

        c = c + 1;
        ifound_types = 1;
    end
    
    answer = strfind(line,'APPROX POSITION XYZ');
    if ~isempty(answer)
        XYZ = sscanf(line(1:61),'%f');
        pos_M = XYZ(1:3);
    end
    
    answer = strfind(line,'ANTENNA: DELTA H/E/N');
    if ~isempty(answer)
        dENU =  sscanf(line(1:61),'%f');
        antoff = [dENU(1:3)];
    end

    answer = strfind(line,'ANT # / TYPE');
    if ~isempty(answer)
        antmod = sscanf(line(21:35),'%c');
        radtyp = sscanf(line(36:40),'%c');
        if (isempty(find(radtyp ~= ' ', 1)))
            radtyp = ' NONE';
        end
        antmod = [antmod radtyp];
    end
    
    answer = strfind(line,'REC # / TYPE / VERS');
    if ~isempty(answer)
        rxmod = sscanf(line(21:40),'%c');
    end
    
    answer = strfind(line,'INTERVAL');
    if ~isempty(answer)
        interval = sscanf(line(1:10),'%f');
    end
    
    %parse next line
    line = fgetl(file);
end

% %apply the antenna offset from the marker (if available)
% if (any(pos_M) && any(antoff))
%     pos_M = local2globalPos(antoff, pos_M);
% end
