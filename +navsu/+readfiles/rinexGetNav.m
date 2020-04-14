function [Eph, iono,glut,leapSeconds] = rinexGetNav(file_nav, constellations)

% SYNTAX:
%   [Eph, iono] = RINEX_get_nav(file_nav, constellations);
%
% INPUT:
%   file_nav = RINEX navigation file
%   constellations = struct with multi-constellation settings (see goGNSS.initConstellation)
%
% OUTPUT:
%   Eph = matrix containing 33 navigation parameters for each satellite
%   iono = matrix containing ionosphere parameters
%
% DESCRIPTION:
%   Parse a RINEX navigation file.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
% Portions of code contributed by Damiano Triglione (2012)
%
% Partially based on RINEXE.M (EASY suite) by Kai Borre
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

% ioparam = 0;
Eph = [];
iono = zeros(8,1);
glut = zeros(4,1);
leapSeconds = nan;

if (nargin < 2 || isempty(constellations)) %then use only GPS as default
    constellations.GPS = struct('numSat', 32, 'enabled', 1, 'indexes', [1:32], 'PRN', [1:32]);
    constellations.nEnabledSat = 32;
    constellations.indexes = constellations.GPS.indexes;
    constellations.PRN     = constellations.GPS.PRN;
end

% fprintf('%s',['Reading RINEX file ' file_nav ': ... ']);

%open navigation file
fid = fopen(file_nav,'rt');

suglFlag = ~isempty(strfind(file_nav,'sugl'));

%read the header
header_end = [];
while (isempty(header_end)) && ~feof(fid)
    %read the line and search the ionosphere labels
    lin = fgetl(fid);
    
    vers_found =  ~isempty(strfind(lin,'RINEX VERSION / TYPE'));
    iono_found = (~isempty(strfind(lin,'ION ALPHA'))) || (~isempty(strfind(lin,'IONOSPHERIC CORR')) && ~isempty(strfind(lin,'BDSA')));
    glut_found = (~isempty(strfind(lin,'GLUT'))) ;
    leap_found = (~isempty(strfind(lin,'LEAP SECONDS'))) ;
    
    %if the ionosphere parameters label was found
    if (vers_found)
        version = str2num(lin(1:9));
    end
    
    if glut_found
        data = textscan(lin(5:end),'%f%f%f%f%*[^\n]');
        glut = cell2mat(data);
    end
    
    if leap_found
        data = textscan(lin(1:end),'%f*[^\n]');
        leapSeconds = cell2mat(data);
    end
    
    %if the ionosphere parameters label was found
    if (iono_found)
        %change flag
        %         ioparam = 1;
        %save the 8 ionosphere parameters
        data = textscan(lin(5:end),'%f%f%f%f%*[^\n]');
        if ~isempty(data{4})
            iono(1) = data{1};
            iono(2) = data{2};
            iono(3) = data{3};
            iono(4) = data{4};
            lin = [];
            while isempty(lin)
                lin = fgetl(fid);
            end
            data = textscan(lin(5:end),'%f%f%f%f%*[^\n]');
            if ~isempty(data{4})
                iono(5) = data{1};
                iono(6) = data{2};
                iono(7) = data{3};
                iono(8) = data{4};
            else
                iono = zeros(8,1);
            end
        end
    end
    
    if glut_found
        
        
    end
    
    header_end = strfind(lin,'END OF HEADER');
end

% GPS TEXTSCAN
formatGPS = ['%1s %d %d %d %d %d %d %d %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n'];

formatSBAS = ['%1s %d %d %d %d %d %d %d %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f'];

formatGLO  = ['%1s %d %d %d %d %d %d %d %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f'];

formatGAL = ['%1s %d %d %d %d %d %d %d %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n'];
formatBDS = ['%1s %d %d %d %d %d %d %d %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f'];

formatQZSS = ['%1s %d %d %d %d %d %d %d %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f'];

formatIRNSS = ['%1s %d %d %d %d %d %d %d %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f %f %f %f \n' ...
    '%f'];
% Look for instances of each constellation entry
fStart = ftell(fid);
allData = textscan(fid,'%s','delimiter','@');
allData = allData{1};

if any(cellfun(@isempty,allData))
    % there's no data here!
    warning(['Error reading ' file_nav])
    fclose(fid);
    return
end

firstChars = cellfun(@(x) x(1), allData, 'un', 0);
constLetters = {'G' 'R' 'E' 'C' 'J' 'I' 'S' 'R2'};
[~,entries] = ismember(firstChars,constLetters);
indEntries = find(entries);
entries(entries == 0) = [];
prnEntries = str2double(cellfun(@(x) x(2:4), allData(indEntries),'un',0));
entries(entries == 2 & prnEntries > 100) = 8; % newest occasional GLOANSS PRN 136 shows up and needs GPS format.

if isempty(entries)
    % there's no data here!
    fclose(fid);
    return
end

% Build a table of sets of outputs by constellation
constTable = [];
constTable = [entries(1) 1 nan];
while true
    % Where is the end of the current section?
    nonSection = find(entries ~= constTable(end,1)) ;
    minSection = min(nonSection(nonSection > constTable(end,2)));
    
    if ~isempty(minSection)
        constTable(end,3) = minSection-1;
        constTable = [constTable; entries(minSection) minSection nan];
    else % we reached the end.
        constTable(end,3) = length(entries);
        break
    end
end


fseek(fid,fStart,'bof');

nConstSections = size(constTable,1);
for idx = 1:nConstSections
    if idx == 14
        'breakpoint';
    end
    const = constTable(idx,1);
    nEntries = constTable(idx,3)-constTable(idx,2)+1;
    
    ephInds = size(Eph,2)+1:(size(Eph,2)+nEntries);
    
    switch const
        case {1 5 6 8} % GPS
            data = textscan(fid,formatGPS,nEntries);
            
            % Only parse and save if constellation is active
            switch constLetters{const}
                case 'G'
                    if (~constellations.GPS.enabled)
                        %                         str = textscan(fid,'%f %f',1,'headerlines',6);
                        continue;
                    end
                case 'E'
                    if (~constellations.Galileo.enabled), continue, end
                case 'C'
                    if (~constellations.BeiDou.enabled), continue, end
                case 'J'
                    if (~constellations.QZSS.enabled), continue, end
                case 'R2'
                    if (~constellations.GLONASS.enabled), continue, end
                case 'I'
                    if (~0), continue, end % IRNSS support ever?
            end
            
            % check if the length of the data is correct
            if sum(cellfun(@length,data)) ~= length(ephInds)*length(data)
                Eph = [];
                iono = [];
                fclose(fid);
                
                fprintf(2,['Error occurred reading ' file_nav ])
                return
            end
            
            svprn = data{2};
            year   = data{3};
            year(year < 20) = year(year < 20)+2000;
            
            month  = data{4};
            day    = data{5};
            hour   = data{6};
            minute = data{7};
            second = data{8};
            
            af0 = data{9};
            af1 = data{10};
            af2 = data{11};
            
            IODE   = data{12};
            crs    = data{13};
            deltan = data{14};
            M0     = data{15};
            
            cuc   = data{16};
            ecc   = data{17};
            cus   = data{18};
            roota = data{19};
            
            toe    = data{20};
            cic    = data{21};
            Omega0 = data{22};
            cis    = data{23};
            
            i0       = data{24};
            crc      = data{25};
            omega    = data{26};
            Omegadot = data{27};
            
            idot       = data{28};
            code_on_L2 = data{29};
            weekno     = data{30};
            L2flag     = data{31};
            
            svaccur  = data{32};
            svhealth = data{33};
            tgd      = data{34};
            iodc     = data{35};
            
            tom = data{36};
            % If time of message is invalid (garbage data input) use 2 hours
            % before ephemeris time.
            if abs(tom) > 86400*7
                tom = toe-7200;
            end
            
            fit_int = data{37};
            
            Eph(1,ephInds)  = svprn;
            Eph(2,ephInds)  = year-2000;
            Eph(3,ephInds)  = month;
            Eph(4,ephInds)  = day;
            Eph(5,ephInds)  = hour;
            Eph(6,ephInds)  = minute;
            Eph(7,ephInds)  = second;
            Eph(8,ephInds)  = af0;
            Eph(9,ephInds)  = af1;
            Eph(10,ephInds) = af2;
            Eph(11,ephInds) = IODE;
            Eph(12,ephInds) = crs;
            Eph(13,ephInds) = deltan;
            Eph(14,ephInds) = M0;
            Eph(15,ephInds) = cuc;
            Eph(16,ephInds) = ecc;
            Eph(17,ephInds) = cus;
            Eph(18,ephInds) = roota;
            Eph(19,ephInds) = toe;
            Eph(20,ephInds) = cic;
            Eph(21,ephInds) = Omega0;
            Eph(22,ephInds) = cis;
            Eph(23,ephInds) = i0;
            Eph(24,ephInds) = crc;
            Eph(25,ephInds) = omega;
            Eph(26,ephInds) = Omegadot;
            Eph(27,ephInds) = idot;
            Eph(28,ephInds) = code_on_L2;
            Eph(29,ephInds) = weekno; if (const == 3); weekno( weekno> 2500) = weekno( weekno> 2500) - 1024; end;
            Eph(30,ephInds) = L2flag;
            Eph(31,ephInds) = svaccur;
            Eph(32,ephInds) = svhealth;
            Eph(33,ephInds) = tgd;
            Eph(34,ephInds) = iodc;
            Eph(35,ephInds) = tom;
            Eph(36,ephInds) = fit_int;
            
        case 2 % GLO
            data = textscan(fid,formatGLO,nEntries);
            
            if (~constellations.GLONASS.enabled), continue, end            
            
            % Parse and add to eph matrix... which needs to be 
            for jdx = 2:length(data)
               Eph(jdx-1,ephInds) = data{jdx}; 
            end
            
        case 3 % GAL
            data = textscan(fid,formatGPS,nEntries);
            
            % Only parse and save if constellation is active
            if (~constellations.Galileo.enabled), continue, end
            
            % check if the length of the data is correct
            if sum(cellfun(@length,data)) ~= length(ephInds)*length(data)
                Eph = [];
                iono = [];
                fclose(fid);
                
                fprintf(2,['Error occurred reading ' file_nav ])
                return
            end
            
            svprn = data{2};
            year   = data{3};
            year(year < 20) = year(year < 20)+2000;
            
            month  = data{4};
            day    = data{5};
            hour   = data{6};
            minute = data{7};
            second = data{8};
            
            af0 = data{9};
            af1 = data{10};
            af2 = data{11};
            
            IODE   = data{12};
            crs    = data{13};
            deltan = data{14};
            M0     = data{15};
            
            cuc   = data{16};
            ecc   = data{17};
            cus   = data{18};
            roota = data{19};
            
            toe    = data{20};
            cic    = data{21};
            Omega0 = data{22};
            cis    = data{23};
            
            i0       = data{24};
            crc      = data{25};
            omega    = data{26};
            Omegadot = data{27};
            
            idot       = data{28};
            code_on_L2 = data{29};
            weekno     = data{30};
            L2flag     = data{31};
            
            svaccur  = data{32};
            svhealth = data{33};
            tgd      = data{34};
            tgd2     = data{35};
            
            tom = data{36};
            % If time of message is invalid (garbage data input) use 2 hours
            % before ephemeris time.
            if abs(tom) > 86400*7
                tom = toe-7200;
            end
            
            fit_int = data{37};
            
            Eph(1,ephInds)  = svprn;
            Eph(2,ephInds)  = year-2000;
            Eph(3,ephInds)  = month;
            Eph(4,ephInds)  = day;
            Eph(5,ephInds)  = hour;
            Eph(6,ephInds)  = minute;
            Eph(7,ephInds)  = second;
            Eph(8,ephInds)  = af0;
            Eph(9,ephInds)  = af1;
            Eph(10,ephInds) = af2;
            Eph(11,ephInds) = IODE;
            Eph(12,ephInds) = crs;
            Eph(13,ephInds) = deltan;
            Eph(14,ephInds) = M0;
            Eph(15,ephInds) = cuc;
            Eph(16,ephInds) = ecc;
            Eph(17,ephInds) = cus;
            Eph(18,ephInds) = roota;
            Eph(19,ephInds) = toe;
            Eph(20,ephInds) = cic;
            Eph(21,ephInds) = Omega0;
            Eph(22,ephInds) = cis;
            Eph(23,ephInds) = i0;
            Eph(24,ephInds) = crc;
            Eph(25,ephInds) = omega;
            Eph(26,ephInds) = Omegadot;
            Eph(27,ephInds) = idot;
            Eph(28,ephInds) = code_on_L2;
            Eph(29,ephInds) = weekno; if (const == 3); weekno( weekno> 2500) = weekno( weekno> 2500) - 1024; end;
            Eph(30,ephInds) = L2flag;
            Eph(31,ephInds) = svaccur;
            Eph(32,ephInds) = svhealth;
            Eph(33,ephInds) = tgd;
            Eph(34,ephInds) = tgd2;
            Eph(35,ephInds) = tom;
            Eph(36,ephInds) = fit_int;
            
            if suglFlag
                suglInfo1 = data{38};
                suglInfo2 = data{39};
                
                Eph(37,ephInds) = suglInfo1;
                Eph(38,ephInds) = suglInfo2;
            end
            
            if 1; %LSB_RECOVERY
                A = Eph';
                piGPS = 3.1415926535898;
                indices = [8:27 33 34];
                %                 scaleFactors = pow2([1 1 1 1 1 piGPS piGPS 1 1 1 1 1 1 piGPS 1 piGPS 1 piGPS piGPS piGPS 1 1], ...
                %                     [-31 -43 -55 0 -5 -43 -31 -29 -33 -29 -19 4 -29 -31 -29 -31 -5 -31 -43 -43 0 0]);
                scaleFactors = pow2([1 1 1 1 1 piGPS piGPS 1 1 1 1 1 1 piGPS 1 piGPS 1 piGPS piGPS piGPS 1 1], ...
                  [-34 -46 -59   0  -6 -43 -31 -31 -33 -31 -19   3 -31 -31 -31 -31  -6 -31 -43 -43  -32 -32]);
 %                [-31 -43 -55   0  -5 -43 -31 -29 -33 -29 -19   4 -29 -31 -29 -31  -5 -31 -43 -43  -31]);
                 % af0 af1 af2 IODE crs dn  M0 cuc   e cus  ra toe cic  O0 cis   i crc   w  Wd  id tgd tgd2
                
                
%                 scaleFactors(end-1:end) = 1e-9*0.1;
                limits = pow2([21 15 7 8 15 15 31 15 32 15 32 16 15 31 15 31 15 31 23 13 10 10]);
                A1 = bsxfun(@rdivide, A(:, indices), scaleFactors);
                A2 = round(A1);
                stderr = nanstd(A2 - A1);
                if max(stderr) > 0.1
%                                         fprintf(2, 'Imprecise %s: %g\n', file_nav, max(stderr));
                end
                %                 idx = any(bsxfun(@ge, abs(A2), limits), 2);
                if any(idx)
                    %                     fprintf(2, 'Over limit %s: %d/%d\n', filename, nnz(idx), length(idx));
                end
                A(:, indices) = bsxfun(@times, A2, scaleFactors);
                Eph = A';
            end
            
        case 4 % BDS
            data = textscan(fid,formatGPS,nEntries);
            
            % Only parse and save if constellation is active
            if (~constellations.BeiDou.enabled), continue, end
            
            % check if the length of the data is correct
            if sum(cellfun(@length,data)) ~= length(ephInds)*length(data)
                Eph = [];
                iono = [];
                fclose(fid);
                
                fprintf(2,['Error occurred reading ' file_nav ])
                return
            end
            
            svprn = data{2};
            year   = data{3};
            year(year < 20) = year(year < 20)+2000;
            
            month  = data{4};
            day    = data{5};
            hour   = data{6};
            minute = data{7};
            second = data{8};
            
            af0 = data{9};
            af1 = data{10};
            af2 = data{11};
            
            IODE   = data{12}; %AODE
            crs    = data{13};
            deltan = data{14};
            M0     = data{15};
            
            cuc   = data{16};
            ecc   = data{17};
            cus   = data{18};
            roota = data{19};
            
            toe    = data{20};
            cic    = data{21};
            Omega0 = data{22};
            cis    = data{23};
            
            i0       = data{24};
            crc      = data{25};
            omega    = data{26};
            Omegadot = data{27};
            
            idot       = data{28};
            code_on_L2 = data{29};
            weekno     = data{30};
            L2flag     = data{31};
            
            svaccur  = data{32};
            svhealth = data{33};
            tgd      = data{34};
            tgd2     = data{35};
            
            tom = data{36};
            %             aoc = data{37};
            % If time of message is invalid (garbage data input) use 2 hours
            % before ephemeris time.
            if abs(tom) > 86400*7
                tom = toe-7200;
            end
            
            fit_int = data{37};
            
            Eph(1,ephInds)  = svprn;
            Eph(2,ephInds)  = year-2000;
            Eph(3,ephInds)  = month;
            Eph(4,ephInds)  = day;
            Eph(5,ephInds)  = hour;
            Eph(6,ephInds)  = minute;
            Eph(7,ephInds)  = second;
            Eph(8,ephInds)  = af0;
            Eph(9,ephInds)  = af1;
            Eph(10,ephInds) = af2;
            Eph(11,ephInds) = IODE;
            Eph(12,ephInds) = crs;
            Eph(13,ephInds) = deltan;
            Eph(14,ephInds) = M0;
            Eph(15,ephInds) = cuc;
            Eph(16,ephInds) = ecc;
            Eph(17,ephInds) = cus;
            Eph(18,ephInds) = roota;
            Eph(19,ephInds) = toe;
            Eph(20,ephInds) = cic;
            Eph(21,ephInds) = Omega0;
            Eph(22,ephInds) = cis;
            Eph(23,ephInds) = i0;
            Eph(24,ephInds) = crc;
            Eph(25,ephInds) = omega;
            Eph(26,ephInds) = Omegadot;
            Eph(27,ephInds) = idot;
            Eph(28,ephInds) = code_on_L2;
            Eph(29,ephInds) = weekno; if (const == 3); weekno( weekno> 2500) = weekno( weekno> 2500) - 1024; end;
            Eph(30,ephInds) = L2flag;
            Eph(31,ephInds) = svaccur;
            Eph(32,ephInds) = svhealth;
            Eph(33,ephInds) = tgd;
            
            % need to convert many tgd2 values
            convInds = find(tgd2 >= 1 & tgd2 <= 1024 & floor(tgd2) == tgd2);
            if ~isempty(convInds)
                tgd2(convInds) = twosComp2dec(dec2bin(tgd2(convInds),10))*1e-9*0.1;
            end

            backwardsScalingInds = find(abs(tgd2)> 1e8);
            if ~isempty(backwardsScalingInds)
                tgd2(backwardsScalingInds) = tgd2(backwardsScalingInds)*(1e-9*0.1)^2;
            end
            
            
            Eph(34,ephInds) = tgd2;
            Eph(35,ephInds) = tom;
            Eph(36,ephInds) = fit_int;
            
            if suglFlag
                suglInfo1 = data{38};
                suglInfo2 = data{39};
                
                Eph(37,ephInds) = suglInfo1;
                Eph(38,ephInds) = suglInfo2;
            end
            
            if 1;%LSB_RECOVERY
                A = Eph';
                piGPS = 3.1415926535898;
                indices = [8:27 33 34];
                %                 scaleFactors = pow2([1 1 1 1 1 piGPS piGPS 1 1 1 1 1 1 piGPS 1 piGPS 1 piGPS piGPS piGPS 1 1], ...
                %                     [-31 -43 -55 0 -5 -43 -31 -29 -33 -29 -19 4 -29 -31 -29 -31 -5 -31 -43 -43 0 0]);
                scaleFactors = pow2([1 1 1 1 1 piGPS piGPS 1 1 1 1 1 1 piGPS 1 piGPS 1 piGPS piGPS piGPS 1 1], ...
                    [-33 -50 -66   0  -6 -43 -31 -31 -33 -31 -19   3 -31 -31 -31 -31  -6 -31 -43 -43   0   0]);
                % af0 af1 af2 IODE crs dn  M0 cuc   e cus  ra toe cic  O0 cis   i crc   w  Wd  id tgd tgd2
                
                
                scaleFactors(end-1:end) = 1e-9*0.1;
                limits = pow2([21 15 7 8 15 15 31 15 32 15 32 16 15 31 15 31 15 31 23 13 10 10]);
                A1 = bsxfun(@rdivide, A(:, indices), scaleFactors);
                A2 = round(A1);
                stderr = nanstd(A2 - A1);
                if max(stderr) > 0.1
%                                         fprintf(2, 'Imprecise %s: %g\n', filename, max(stderr));
                end
                %                 idx = any(bsxfun(@ge, abs(A2), limits), 2);
                if any(idx)
                    %                     fprintf(2, 'Over limit %s: %d/%d\n', filename, nnz(idx), length(idx));
                end
                A(:, indices) = bsxfun(@times, A2, scaleFactors);
                Eph = A';
                %                 A(idx, :) = [];
            end
            
            
        case 7 % SBAS
            data = textscan(fid,formatSBAS,nEntries);
            if (~constellations.SBAS.enabled), continue, end
    end
    
end

fclose(fid);
