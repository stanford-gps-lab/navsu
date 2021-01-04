function [Eph,ionoCorr,timeCorr,leapSecs] = rinexGetNav(file_nav, constellations)
% SYNTAX:
%   [Eph, ionoCorr, timeCorr, leapSecs] = RINEX_get_nav(file_nav, constellations)
%
% INPUT:
%   file_nav = RINEX navigation file
%   constellations = struct with multi-constellation settings (see goGNSS.initConstellation)
%
% OUTPUT:
%   Eph = matrix containing 33 navigation parameters for each satellite
%   ionoCorr = 8x1 matrix containing iono params (ai0-ai2 plus a blank [GAL];
%              or alpha0--alpha3 and/or beta0--beta3 [GPS/QZS/BDS/IRNSS]) [OPTIONAL]
%              Note: time mark and ID of satellite  [RINEX 3 only] [OPTIONAL]
%   timeCorr = matrix containing time system correction parameters [OPTIONAL]
%   leapSeconds = current #of leap seconds in UTC relative to GNSS system time [OPTIONAL]
%
%   RINEX specification allows params marked [OPTIONAL] to be omitted from
%   header. In such cases, zeros are returned in these outputs.
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
%    This program is free software: you can redistribute it and/or modify`
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

%% General setup

% verbose output -- flag should be assigned in calling function
global DEBUG;

Eph = [];
[vers_found, timeCorrR2_found, timeCorrR3_found] = deal(0);
[timeCorr, ionoCorr] = deal(NaN);
ionoCorrCoeffs = NaN(8,1);
leapSecs = NaN(5,1);

if (nargin < 2 || isempty(constellations)) %then use only GPS as default
    constellations.GPS = struct('numSat', 32, 'enabled', 1, 'indexes', 1:32, 'PRN', 1:32);
    constellations.nEnabledSat = 32;
    constellations.indexes = constellations.GPS.indexes;
    constellations.PRN     = constellations.GPS.PRN;
end


%% Read contents of navigation file

if DEBUG,
    fprintf('Reading RINEX file %s ...\n\n', file_nav);
end

fid = fopen(file_nav,'rt');
% Note: by default, TEXTSCAN strips leading whitespace. Since many RINEX
%       format specifiers start with mandatory skips (e.g. 3X), force it
%       to preserve leading whitespace when making this call. This means
%       that all string parsing format specifiers throughout this script
%       must precisely match those in the RINEX spec, or be tolerant to
%       subsequent TEXTSCAN calls below.
allData = textscan(fid,'%s','whitespace','','Delimiter','\n');
allData = allData{1};
fclose(fid);

%perform some simple checks
if any(cellfun(@isempty,allData))
    warning('Error reading %s -- unexpected blank line(s) found!\n', file_nav);
    return
else
    header_end = find(cellfun(@(x) ~isempty(x), strfind(allData,'END OF HEADER')));
    if isempty(header_end) 
        warning('Error reading %s -- end of header not found!\n', file_nav);
        return
    end
    
    if length(header_end) > 1
        warning('Error reading %s -- multiple headers found!\n', file_nav);
       return 
    end
end

suglFlag = ~isempty(strfind(file_nav,'sugl'));


%% Parse the header

for idx = 1:header_end
    
    lin = allData{idx};
    
    %%% RINEX version and file type
    if ~isempty(strfind(lin, 'RINEX VERSION / TYPE'))
        vers_found = 1; 
        rinexVersion = cell2mat(textscan(lin(1:9), '%9.2f')); 
        rinexType = lin(21);
        
        % sanity check version and system type
        if floor(rinexVersion) == 2, % Known versions are 2.00 -- 2.11 (May 2020)
            if ~ismember(rinexType, 'NGH'), % N=GPS[Table A3]; G=GLONASS[Table A10]; H=SBAS/Geo [Table A15]
                error(['Error in header of %s -- RINEX %.2f file type = %s is undefined ' ...
                    '(allowed: N [GPS], G [GLONASS], H [SBAS/Geo]).\n'], ...
                    file_nav, rinexVersion, rinexType);
            end
        elseif floor(rinexVersion) == 3, % Known versions are 3.00 -- 3.04 (May 2020)
            if rinexType ~= 'N',
                error(['Error in header of %s -- RINEX %.2f File Type = %s is undefined ' ...
                    '(should be ''N'' for navigation message file).\n'], ...
                    file_nav, rinexVersion, rinexType);
            end
            rinexSatelliteSystem = lin(41); % single char; one of G/R/E/J/C/I/S/M
            if ~ismember(rinexSatelliteSystem, 'GREJCISM'),
                error(['Error in header of %s -- RINEX %.2f satellite system code = %s is undefined ' ...
                    '(should be one of G / R / E / J / C / I / S / M)'], ...
                    file_nav, rinexVersion, rinexSatelliteSystem);
            end
        else
            error('Error in header of %s -- RINEX version %.2f unknown!', file_nav, rinexVersion);
        end
        
    end
    
    %%% Time system corrections [OPTIONAL]
    %
    % RINEX 2: (Table A3) Time corrections must be marked with the header
    %          label (columns 61-80) 'DELTA-UTC: A0,A1,T,W'. Note that GLUT
    %          is NOT defined in RINEX 2.
    if ~isempty(strfind(lin,'DELTA-UTC: A0,A1,T,W')),
        timeCorrR2_found = 1;
        data = textscan(lin, '%19.12f %19.12f %9d %9d %*[^\n]');
        timeCorrType = 'DUTC'; % analogous to GPUT in RINEX 3.x?
        [a0, a1, T, W] = data{:};
        [timeSource, utcIdentifier] = deal(NaN); % not specified in RINEX 2
    end
    % RINEX 3: (Table A5) Time corrections must be marked with the header
    %          label (columns 61-80) 'TIME SYSTEM CORR'. Although the spec
    %          does not say so explicitly (unlike iono corrections, which
    %          it notes *can* be provided in multiple messages), it appears
    %          that multiple time correction lines may legally appear in the
    %          header. Each such line indicates corrections between two
    %          specific time systems, taken pairwise (UTC, GPS, GAL, GLO,
    %          QZSS, IRNSS, and SBAS). The first field of each such line
    %          must be marked with a 4-character identifier specifying the 
    %          two systems related by the parameters given on that line. All
    %          such lines contain the same four parameters: a0, a1, T and W.
    %
    %          (Tables A7, A9, A11, A15) In the example TIME SYSTEM CORR
    %          header lines, the 5th and 6th fields ("S", a 5-char service
    %          provider [for SBAS, derived from message type 17 or 12;
    %          otherwise, the ID of the SV broadcasting the message]; and
    %          "U", a 2-digit-integer UTC identifier, respectively) are
    %          omitted in RINEX 3.3; only the SBAS example (Table A17)
    %          includes these fields. In RINEX 3.4, all examples are revised
    %          to include these fields.
    timeCorrIdx = strfind(lin, 'TIME SYSTEM CORR'); % handle (possibly missing)
                                                    % service provider, UTC
                                                    % ID without confusing
                                                    % TEXTSCAN
    if ~isempty(timeCorrIdx),
        timeCorrR3_found = 1;
        data = textscan(lin(1:timeCorrIdx-1), '%4c%17.10f%16.9f%7d%5d%5c%2d%*[^\n]');
        [timeCorrType, timeSource] = data{[1 6]}; % char/string
        [a0, a1, T, W, utcIdentifier] = data{[2:5 7]}; % float/int
    end
    % Whether RINEX 2 or 3, assemble "Time System Corrections" struct(s)
    if (timeCorrR2_found || timeCorrR3_found),
        dummy = struct( ...
            'timeCorrType', timeCorrType, ...
            'a0', a0, 'a1', a1, 'T', T, 'W', W, ...
            'timeSource', timeSource, 'UTCID', utcIdentifier );
        if ~isstruct(timeCorr),
            timeCorr = dummy; % first line encountered
        else
            timeCorr = [timeCorr dummy]; %#ok<AGROW> % subsequent TIME SYSTEM CORR lines
        end
        [timeCorrR2_found, timeCorrR3_found] = deal(0); % avoid extraneous warnings
    end
    
    
    %%% Iono corrections [OPTIONAL]
    %
    % RINEX 2: (Table A3) A full set of iono correction params requires both
    %          lines marked "ION ALPHA" (for coeffs. alpha0 - alpha3) AND
    %          "ION BETA" (for coeffs. beta0 - beta3). Check for presence of
    %          both in populating iono coefficient vector, otherwise print
    %          warning and fill remainder with zeros.
    if ~isempty(strfind(lin,'ION ALPHA')),
        data = textscan(lin(3:end), '%12.4f%12.4f%12.4f%12.4f%*[^\n]');
        ionoCorrCoeffs(1:4) = cell2mat(data); % A0-A3
        ionoCorrType = 'RINEX2_A0-B3';
        timeMark = ''; ionoSVID = NaN; % not specified in RINEX 2
    end
    if ~isempty(strfind(lin,'ION BETA')),
        data = textscan(lin(3:end), '%12.4f%12.4f%12.4f%12.4f%*[^\n]');
        ionoCorrCoeffs(5:8) = cell2mat(data); % B0-B3
        timeMark = ''; ionoSVID = NaN; % not specified in RINEX 2
    end
    % RINEX 3: (Table A5) Multiple sets of iono correction params can appear
    %          in a single header block. Each such line must be marked with
    %          the header label (col. 61-80) 'IONOSPHERIC CORR'. Two more
    %          values, time mark and SV ID, may appear after the iono
    %          correction parameters on any line, and are mandatory for BDS
    %          (BeiDou) but optional for the other constellations.
    ionoCorrIdx = strfind(lin, 'IONOSPHERIC CORR'); % handle (possibly missing)
                                                    % time mark and UTC ID
                                                    % without confusing TEXTSCAN
    if ~isempty(ionoCorrIdx),
        
        data = textscan(lin(1:ionoCorrIdx-1), '%4s%12.4f%12.4f%12.4f%12.4f%c%2d%*[^\n]');
        
        switch cell2mat(data{1})
            case 'GAL'
                if ~any(isnan(ionoCorrCoeffs)), % no longer NaNs once an IONO CORR line is parsed
                    warning([ ...
                        'In input file %s: \n' ...
                        'Multiple Galileo IONOSPHERIC CORR lines found -- overwriting ' ...
                        'previously parsed parameters. Need to \nmodify output argument '...
                        'struct(s) to handle parameters from multiple systems!\n' ]);
                end
                % 'GAL', ai0 - ai2, then a blank
                ionoCorrType = cell2mat(data{1}); 
                ionoCorrCoeffs(1:3) = deal(cell2mat(data(2:4)));
                ionoCorrCoeffs(4:8) = 0;
                timeMark = data(6); ionoSVID = data{7}; % transmission time, SVID (mandatory for BDS, optional for other systems)                
            case {'GPSA','QZSA','BDSA','IRNA'}
                if ~any(isnan(ionoCorrCoeffs(1:4))), % no longer NaNs once an IONO CORR line is parsed
                    warning(['In input file %s:\n' ...
                        'Multiple GPSA/QZSA/BDSA/IRNA IONOSPHERIC CORR lines found! ' ...
                        'Overwriting previously parsed parameters.\n' ...
                        'Need to modify output argument struct(s) to handle parameters from multiple systems!\n'], ...
                        file_nav);
                end
                % 'xxxA', then alpha0 - alpha3
                ionoCorrType = cell2mat(data{1}); 
                ionoCorrCoeffs(1:4) = deal(cell2mat(data(2:5)));
                % these are required for BDS, optional for other systems:
                if exist('timeMark','var') && ~strcmp(timeMark, data{6})
                    warning(['Iono correction time mark in line with label ''%s'' has different ' ...
                        'time mark than previously parsed value.\n'], ionoCorrType);
                end
                if exist('ionoSVID','var') && ~isempty(ionoSVID) && (ionoSVID~=data{7})
                    warning(['Iono correction time mark in line with label ''%s'' has different ' ...
                        'SVID (%d) than previously parsed value (%d).\n'], ionoCorrType, data{7}, ionoSVID);
                end
                timeMark = data(6); ionoSVID = data{7}; % transmission time, SVID (mandatory for BDS, optional for other systems)                
                if (strcmp(data{1},'BDSA')) && ~all(data{6:7})
                    warning('Error reading %s -- missing mandatory time mark and/or SVID field in BeiDou line!\n', file_nav);
                end
            case {'GPSB','QZSB','BDSB','IRNB'}
                if ~any(isnan(ionoCorrCoeffs(5:8))) % no longer NaNs once an IONO CORR line is parsed
                    warning(['In input file %s:\n' ...
                        'Multiple GPSB/QZSB/BDSB/IRNB IONOSPHERIC CORR lines found! ' ...
                        'Overwriting previously parsed parameters.\n' ...
                        'Need to modify output argument struct(s) to handle parameters from multiple systems!\n'], ...
                        file_nav);
                end
                % 'xxxB', then beta0 - beta3
                ionoCorrType = cell2mat(data{1});
                ionoCorrCoeffs(5:8) = deal(cell2mat(data(2:5)));
                % these are required for BDS, optional for other systems:
                if exist('timeMark','var') && ~strcmp(timeMark,data{6})
                    warning(['Iono correction time mark in line with label ''%s'' has different ' ...
                        'time mark than previously parsed value.\n'], ionoCorrType);
                end
                if exist('ionoSVID','var') && ~isempty(ionoSVID) && (ionoSVID~=data{7})
                    warning(['Iono correction time mark in line with label ''%s'' has different ' ...
                        'SVID (%d) than previously parsed value (%d).\n'], ionoCorrType, data{7}, ionoSVID);
                end
                timeMark = data(6); ionoSVID = data{7}; % transmission time, SVID (mandatory for BDS, optional for other systems)                
                if (strcmp(data{1},'BDSB')) && ~all(data{6:7})
                    warning('Error reading %s -- missing mandatory time mark and/or SVID field in BeiDou line!\n', file_nav);
                end
            otherwise
                warning('Error reading %s -- unexpected field(s) in IONO line:\n   %s\n', file_nav, lin);
                return;
        end % switch

    end % if ~isempty(ionoCorrIdx
    
    % Whether RINEX 2 or 3, assemble "Iono Corrections" struct
    if all(~isnan(ionoCorrCoeffs))
        if ~strcmp(ionoCorrType, 'RINEX2_A0-B3')
            ionoCorrType = ionoCorrType(1:3); % strip off 'A'/'B'
        end
        dummy = struct( ...
            'ionoCorrType', ionoCorrType, ... 
            'ionoCorrCoeffs', ionoCorrCoeffs, ...
            'timeMark', timeMark, 'ionoSVID', ionoSVID);
        ionoCorrCoeffs = nan(8,1); % avoid duplicating parsed data into iono struct next time around
        if ~isstruct(ionoCorr)
            ionoCorr = dummy; % first line encountered
            clear dummy;
        else
            ionoCorr = [ionoCorr dummy]; %#ok<AGROW> % subsequent IONOSPHERIC CORR lines
            clear dummy;
        end
    end
    
    %%% Leap second [OPTIONAL]
    leapInd = strfind(lin,'LEAP SECONDS');
    if ~isempty(leapInd)
        leap_found = 1; %#ok<NASGU>
        data = textscan(lin(1:leapInd-1),'%6d%6d%6d%6d%3c%*[^\n]');
        dataIdx = find(~cellfun(@isempty,data(1:end-1))); %if non-empty, last element (time sys ID) is a string
        leapSecs(dataIdx) = cell2mat(data(dataIdx));
        if ~cellfun(@isempty, data(end))
            % No obvious place to return time system ID string found in
            % LEAP SECOND header line, so at just print warning for now.
            warning('Discarding time system identifier ''%s'', found in LEAP SECOND header line.\n', data{end});
        end
    end
    
end % for idx = 1:header_end

% Enforce mandatory header version info (but ignore contents of PGM / RUN BY / DATE
% for now, since that info is not needed for our GNSS calculations)
if ~vers_found
    error('Error in header of %s -- no RINEX version number found.\n', file_nav);
end
    
% Summarize parsed header
if DEBUG
    fprintf('Parsed header:\n\n');
    fprintf('Version = %.2f\nFile type = ''%s''\nSatellite system = ''%s''\n', ...
        rinexVersion, rinexType, rinexSatelliteSystem);
    fprintf('\ntimeCorr =\n\n');
    if isstruct(timeCorr)
        % use 'AsArray' because some files cause this struct to have fields with different numbers of rows
        disp(struct2table(timeCorr, 'AsArray', true));
    else
        timeCorr, %#ok<NOPRT>
    end
    fprintf('\nionoCorr =\n\n');
    if isstruct(ionoCorr)
        % use 'AsArray' because some files cause this struct to have fields with different numbers of rows
        disp(struct2table(ionoCorr, 'AsArray', true));
        fprintf('ionoCorrCoeffs =\n\n%s\n', ...
            evalc('disp(cell2mat(struct2table(ionoCorr,''AsArray'',true).ionoCorrCoeffs''))'));
    else
        ionoCorr, %#ok<NOPRT>
    end
    fprintf('\nleapSecs = %s\n', mat2str(leapSecs'));
end % if DEBUG


%% Set up TEXTSCAN format specifiers for defined GNSS entries

constData = struct( ...
    'constLetter', {'G'  'E'  'R'  'J'  'C'  'S'  'I'  'R2'}, ...
    'blockLines', {8 8 4 8 8 4 8 4}, ...
    'formatString', cell(1,8), ... % placeholder; filled in below
    'enabled', cell(1,8) ... % placeholder; allows programmatic addressing in addition to "constellations.*.enabled" calls
);

recordFormat = '%f ';   % This previously '%19.12f ', but not all RINEX creators adhere to the standard

% GPS Data Record (RINEX 3.03, Table A6) -- System Identifier "G"
constData(1).formatString = [ ...
    '%1s %2.2d %4d ' repmat('%2.2d ',1,4) repmat(recordFormat,1,4) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1; see note "***"
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
    repmat(recordFormat,1,4) ... % broadcast orbit - 4
    repmat(recordFormat,1,4) ... % broadcast orbit - 5
    repmat(recordFormat,1,4) ... % broadcast orbit - 6
    repmat(recordFormat,1,4) ... % broadcast orbit - 7
    ];
constData(1).enabled = constellations.GPS.enabled;

% Galileo Data Record (RINEX 3.03, Table A8) -- System Identifier "E"
constData(2).formatString = [ ...
    '%1s %2.2d %4d ' repmat('%2.2d ',1,5) repmat(recordFormat,1,3) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1 -- see note "***"
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
    repmat(recordFormat,1,4) ... % broadcast orbit - 4
    repmat(recordFormat,1,4) ... % broadcast orbit - 5 -- see note "****"
    repmat(recordFormat,1,4) ... % broadcast orbit - 6 -- see note "*****"
    repmat(recordFormat,1,4) ... % broadcast orbit - 7 -- see note "**"
    ];
constData(2).enabled = constellations.Galileo.enabled;

% GLONASS Data Record (RINEX 3.03, Table A10) -- System Identifier "R"
constData(3).formatString = [ ...
    '%1s %2.2d %4d ' repmat('%2.2d ',1,5) repmat(recordFormat,1,3) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
    ];
constData(3).enabled = constellations.GLONASS.enabled;

% QZSS Data Record (RINEX 3.03, Table A12) -- System Identifier "J"
constData(4).formatString = [ ...
    '%1s %2d %4d ' repmat('%2d ',1,5) repmat(recordFormat,1,3) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
    repmat(recordFormat,1,4) ... % broadcast orbit - 4
    repmat(recordFormat,1,4) ... % broadcast orbit - 5
    repmat(recordFormat,1,4) ... % broadcast orbit - 6
    repmat(recordFormat,1,4) ... % broadcast orbit - 7 -- see note "**"
    ];
constData(4).enabled = constellations.QZSS.enabled;

% BDS (BeiDou) Data Record (RINEX 3.03, Table A14) -- System Identifier "C"
constData(5).formatString = [ ...
    '%1s %2.2d %4d ' repmat('%2.2d ',1,5) repmat(recordFormat,1,3) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1 -- see note "**"
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
    repmat(recordFormat,1,4) ... % broadcast orbit - 4
    repmat(recordFormat,1,4) ... % broadcast orbit - 5 -- see note "***"
    repmat(recordFormat,1,4) ... % broadcast orbit - 6
    repmat(recordFormat,1,4) ... % broadcast orbit - 7 -- see note "****"
    ];
constData(5).enabled = constellations.BeiDou.enabled;

% SBAS Data Record (RINEX 3.03, Table A16) -- System Identifier "S"
constData(6).formatString = [ ...
    '%1s %2.2d %4d ' repmat('%2.2d ',1,5) repmat(recordFormat,1,3) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
];
constData(6).enabled = constellations.SBAS.enabled;

% IRNSS Data Record (RINEX 3.03, Table A18) -- System Identifier "I"
constData(7).formatString = [ ...
    '%1s %2.2d %4d ' repmat('%2.2d ',1,5) repmat(recordFormat,1,3) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1 -- see note "***"
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
    repmat(recordFormat,1,4) ... % broadcast orbit - 4
    repmat(recordFormat,1,4) ... % broadcast orbit - 5
    repmat(recordFormat,1,4) ... % broadcast orbit - 6
    repmat(recordFormat,1,4) ... % broadcast orbit - 7 -- see note "**"
    ];
constData(7).enabled = 0; % As of May 2020, constellations.IRNSS.enabled not defined yet...

% GLONASS2 Data Record (RINEX 3.03, Table A10[??]) -- "R2" [??]
constData(8).formatString = [ ...
    '%1s %2.2d %4d ' repmat('%2.2d ',1,5) repmat(recordFormat,1,3) ... % SV / Epoch / SV Clk -- see note "*"
    repmat(recordFormat,1,4) ... % broadcast orbit - 1
    repmat(recordFormat,1,4) ... % broadcast orbit - 2
    repmat(recordFormat,1,4) ... % broadcast orbit - 3
    ];
constData(8).enabled = constellations.GLONASS.enabled;


%% Parse navigation data lines

% Look for instances of each constellation entry, skipping the header text
firstChars = cellfun(@(x) x(1), allData((header_end+1):end), 'un', 0);

% check GNSS identifiers (first character of each data block)
if ~any(cellfun(@(x) ismember(x,[constData(1:end-1).constLetter]), firstChars)) 
    if floor(rinexVersion) == 2  % RINEX 2.x
        if DEBUG, warning('No constellation identifiers in first column of input. (This is LEGAL in RINEX 2.x)'); end
        if rinexType == 'N'      % GPS
            dummy = 'G';
        elseif rinexType == 'G'  % GLONASS
            dummy = 'R';
        else
            dummy = 'S';          % SBAS/GEO
        end
    else % RINEX v3.x
        if rinexType == 'M'
            error(['Error in header of %s -- File Type = M (mixed constellation data),\n' ...
                'so no way to tell which GNSS each parameter block refers to.\n\n'], file_nav);
        else
            dummy = rinexType; % assume all parameter blocks are for the indicated GNSS
        end            
    end
    % Insert appropriate identifiers, padding first column of remaining lines w/ spaces
    if DEBUG, warning('Inserting ''%c'' at start of each parameter block.\n\n', dummy); end
    firstLines = find(cellfun(@(x) ~isspace(x(2)), allData((header_end+1):end))); % first line of each parameter block; insert identifier here
    firstChars = repmat({' '}, size(allData,1)-header_end, 1); % one-space column for padding
    firstChars(firstLines) = repmat({dummy}, size(firstLines,1), 1); 
    allData((header_end+1):end) = strcat(firstChars, allData((header_end+1):end));
end
        
%constLetters = {'G' 'E' 'R' 'J' 'C' 'S' 'I' 'R2'};
[~,entries] = ismember(firstChars,{constData.constLetter});
entries(entries == 0) = [];
indEntries = find(entries);
prnEntries = str2double(cellfun(@(x) x(2:4), allData(indEntries), 'un', 0)); %#ok<FNDSB>
entries(entries == 2 & prnEntries > 100) = 8; % newest occasional GLOANSS PRN 136 shows up and needs GPS format.

% Find number of SVs listed for each constellation
entryCounts = num2cell(hist(entries,1:8)); %#ok<HIST>
%entryCounts = num2cell(histogram(entries,'BinLimits',[1,8],'BinMethod','integers').Values); %recommended for R2014a or newer
[constData.numSats] = entryCounts{:};

% Check if we have any data to work with...
if isempty(entries)
    error('In input file %s: no GNSS navigation data entries found!\n', file_nav);
end

for constIdx = 1:length(constData)
    
    if DEBUG, fprintf('constIdx = %d, constLetter = "%s" ...\n', constIdx, constData(constIdx).constLetter); end %#ok<*UNRCH>
    
    % Skip constellations not requested, or with no entries in input file
    if ~constData(constIdx).numSats
        if DEBUG, fprintf('Skipping constellation %s -- no entries in input file...\n', constData(constIdx).constLetter); end
        continue;
    % % Include the following block for slighly more efficient parsing; the short-circuit logic is
    % % already handled by the "if(~constellations.*.enabled), continue; end" blocks below, but 
    % % this part will avoid repeatedly parsing input lines that would get discarded later anyway
    elseif ~constData(constIdx).enabled
        if DEBUG, fprintf('Skipping constellation %s -- processing not requested by calling function...\n', constData(constIdx).constLetter); end
        continue;
    end
    
    % Row indices of FIRST lines of all SV data blocks for this constellation
    constData(constIdx).firstLinesIdx = header_end + find(strcmp(constData(constIdx).constLetter, firstChars));
    
    % Number of satellites for this constellation
    %%% UNNECESSARY -- done above using hist() for all constellations in a single step
    % constData(constIdx).numSats = length(constData(constIdx).firstLinesIdx);
    
    % Row indices of ALL SV data blocks for this constellation, based upon FIRST line indices found above
    dummy = arrayfun(@(x) x+(0:(constData(constIdx).blockLines - 1)), constData(constIdx).firstLinesIdx, 'UniformOutput', false);
    constData(constIdx).allLines = [dummy{:}];
    
    % Reshape all entries for this constellation into correct format for TEXTSCAN,
    % padding with spaces to handle longer lines (e.g. those containing data in 'spare' fields
    dummy2 = join( reshape( allData(constData(constIdx).allLines), constData(constIdx).blockLines, constData(constIdx).numSats )' );
    maxLen = max(cellfun(@(x) size(x,2), dummy2)); % widest row
    dummy2padded = cellfun(@(x) [x repmat(' ', 1, maxLen-size(x,2))], dummy2, 'UniformOutput', false);
    
%     if length(unique(cellfun(@(x) size(x, 2), dummy2))) > 1,
%         warning('One or more lines in input FAILS because not all rows of dummy2 are of equal lengths... PAUSED');
%         pause;
%     end
    dummy3 = [ vertcat(dummy2padded{:}) repmat(char(13), constData(constIdx).numSats, 1) ];
    
    % Parse all entries for this constellation into cell array using appropriate format string
    data = textscan(reshape(dummy3', 1, []), constData(1).formatString);

    switch constIdx
        
        case {1 4 7 8} % GPS/QZSS/IRNSS/GLONASSv2 (RINEX markers: "G"/"J"/"I"/"R2")
            
            % Only parse and save if constellation is active
            switch constData(constIdx).constLetter
                case 'G'
                    if (~constellations.GPS.enabled), continue, end
                case 'J'
                    if (~constellations.QZSS.enabled), continue, end
                case 'R2'
                    if (~constellations.GLONASS.enabled), continue, end
                case 'I'
                    if (~0), continue, end % IRNSS support ever?
            end
          
            svprn  = data{2};
            year   = data{3};
            year(year<80) = year(year<80)+2000; % See RINEX 3 section 6.10
            year(year<=99) = year(year<=99)+1900; % See RINEX 3 section 6.10
            
            month  = data{4};
            day    = data{5};
            hour   = data{6};
            minute = data{7};
            second = data{8};
            
            af0    = data{9};
            af1    = data{10};
            af2    = data{11};
            
            IODE   = data{12};
            crs    = data{13};
            deltan = data{14};
            M0     = data{15};
            
            cuc    = data{16};
            ecc    = data{17};
            cus    = data{18};
            roota  = data{19};
            
            toe    = data{20};
            cic    = data{21};
            Omega0 = data{22};
            cis    = data{23};
            
            i0         = data{24};
            crc        = data{25};
            omega      = data{26};
            Omegadot   = data{27};
            
            idot       = data{28};
            code_on_L2 = data{29};
            weekno     = data{30};
            L2flag     = data{31};
            
            svaccur    = data{32};
            svhealth   = data{33};
            tgd        = data{34};
            iodc       = data{35};
            
            tom        = data{36};
            % If time of message is invalid (garbage data input) use 2 hours
            % before ephemeris time. *** See also RINEX 3.4, section 6.1.3,
            % which prescribes adjusting this ToM by -604800 relative to the
            % received value, resulting in a negative value, for messages
            % broadcast near the midnight Saturday/Sunday (UTC) boundary,
            % when ToE and ToC already refer to the following week.
            if abs(tom) > 86400*7
                tom = toe-7200;
            end
            
            fit_int    = data{37};
            
            Eph(1,:)  = svprn;
            Eph(2,:)  = year-2000;
            Eph(3,:)  = month;
            Eph(4,:)  = day;
            Eph(5,:)  = hour;
            Eph(6,:)  = minute;
            Eph(7,:)  = second;
            Eph(8,:)  = af0;
            Eph(9,:)  = af1;
            Eph(10,:) = af2;
            Eph(11,:) = IODE;
            Eph(12,:) = crs;
            Eph(13,:) = deltan;
            Eph(14,:) = M0;
            Eph(15,:) = cuc;
            Eph(16,:) = ecc;
            Eph(17,:) = cus;
            Eph(18,:) = roota;
            Eph(19,:) = toe;
            Eph(20,:) = cic;
            Eph(21,:) = Omega0;
            Eph(22,:) = cis;
            Eph(23,:) = i0;
            Eph(24,:) = crc;
            Eph(25,:) = omega;
            Eph(26,:) = Omegadot;
            Eph(27,:) = idot;
            Eph(28,:) = code_on_L2;
            Eph(29,:) = weekno;
            Eph(30,:) = L2flag;
            Eph(31,:) = svaccur;
            Eph(32,:) = svhealth;
            Eph(33,:) = tgd;
            Eph(34,:) = iodc;
            Eph(35,:) = tom;
            Eph(36,:) = fit_int;         
            
            
        case 2 % GAL (RINEX marker: "E")
            
            % Only parse and save if constellation is active
            if (~constellations.Galileo.enabled), continue, end
                       
            svprn  = data{2};
            year   = data{3};
            year(year<80) = year(year<80)+2000; % See RINEX 3 section 6.10
            year(year<=99) = year(year<=99)+1900; % See RINEX 3 section 6.10
            
            month  = data{4};
            day    = data{5};
            hour   = data{6};
            minute = data{7};
            second = data{8};
            
            af0 = data{9};
            af1 = data{10};
            af2 = data{11};
            
            IODnav = data{12}; % Analogous to IODE of other GNSSes
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
            % Handle week number rollovers
            weekno(weekno > 2500) = weekno(weekno > 2500) - 1024;
            
            
            if size(dummy3,2) < 580
                % Some files don't put any data in the spare field (rather
                % than using zeros), which offsets the outputs
               noSpareData = 1;
               L2flag = zeros(size(data{31}));
            else
               noSpareData = 0; 
               L2flag     = data{31};
            end            
            
            svaccur  = data{32-noSpareData};
            svhealth = data{33-noSpareData};
            tgd      = data{34-noSpareData};
            tgd2     = data{35-noSpareData};
            
            tom = data{36-noSpareData};
            % If time of message is invalid (garbage data input) use 2 hours
            % before ephemeris time. *** See also RINEX 3.4, section 8.3.3,
            % which prescribes adjusting ToM by +/-604800 relative to the
            % received value, rather than 7200 sec before ToE, as here.
            if abs(tom) > 86400*7
                tom = toe-7200;
            end
            
            fit_int = data{37};
            

            Eph(1,:)  = svprn;
            Eph(2,:)  = year-2000;
            Eph(3,:)  = month;
            Eph(4,:)  = day;
            Eph(5,:)  = hour;
            Eph(6,:)  = minute;
            Eph(7,:)  = second;
            Eph(8,:)  = af0;
            Eph(9,:)  = af1;
            Eph(10,:) = af2;
            Eph(11,:) = IODnav;
            Eph(12,:) = crs;
            Eph(13,:) = deltan;
            Eph(14,:) = M0;
            Eph(15,:) = cuc;
            Eph(16,:) = ecc;
            Eph(17,:) = cus;
            Eph(18,:) = roota;
            Eph(19,:) = toe;
            Eph(20,:) = cic;
            Eph(21,:) = Omega0;
            Eph(22,:) = cis;
            Eph(23,:) = i0;
            Eph(24,:) = crc;
            Eph(25,:) = omega;
            Eph(26,:) = Omegadot;
            Eph(27,:) = idot;
            Eph(28,:) = code_on_L2;
            Eph(29,:) = weekno;
            Eph(30,:) = L2flag;
            Eph(31,:) = svaccur;
            Eph(32,:) = svhealth;
            Eph(33,:) = tgd;
            Eph(34,:) = tgd2;
            Eph(35,:) = tom;
            Eph(36,:) = fit_int;
                      
            if suglFlag
                suglInfo1 = data{38};
                suglInfo2 = data{39};
                
                Eph(37,:) = suglInfo1;
                Eph(38,:) = suglInfo2;
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
            
            
        case 3 % GLOv1 (RINEX marker: "R")
            
            if (~constellations.GLONASS.enabled), continue, end            
                        
            % Parse and add to eph matrix... which needs to be
            %
            % %XXX needs to be what? the above (incomplete) comment is in
            % %XXX original GitHub checkout as of April 2020
            for jdx = 2:length(data)
               Eph(jdx-1,:) = data{jdx}; 
            end
            
            
        case 5 % BDS (RINEX marker: "C")
            
            % Only parse and save if constellation is active
             if (~constellations.BeiDou.enabled), continue, end
                        
            svprn = data{2};
            year   = data{3};
            year(year<80) = year(year<80)+2000; % See RINEX 3 section 6.10
            year(year<=99) = year(year<=99)+1900; % See RINEX 3 section 6.10
            
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
            
            Eph(1,:)  = svprn;
            Eph(2,:)  = year-2000;
            Eph(3,:)  = month;
            Eph(4,:)  = day;
            Eph(5,:)  = hour;
            Eph(6,:)  = minute;
            Eph(7,:)  = second;
            Eph(8,:)  = af0;
            Eph(9,:)  = af1;
            Eph(10,:) = af2;
            Eph(11,:) = IODE;
            Eph(12,:) = crs;
            Eph(13,:) = deltan;
            Eph(14,:) = M0;
            Eph(15,:) = cuc;
            Eph(16,:) = ecc;
            Eph(17,:) = cus;
            Eph(18,:) = roota;
            Eph(19,:) = toe;
            Eph(20,:) = cic;
            Eph(21,:) = Omega0;
            Eph(22,:) = cis;
            Eph(23,:) = i0;
            Eph(24,:) = crc;
            Eph(25,:) = omega;
            Eph(26,:) = Omegadot;
            Eph(27,:) = idot;
            Eph(28,:) = code_on_L2;
            Eph(29,:) = weekno;
            Eph(30,:) = L2flag;
            Eph(31,:) = svaccur;
            Eph(32,:) = svhealth;
            Eph(33,:) = tgd;
 
            % need to convert many tgd2 values
            convInds = find(tgd2 >= 1 & tgd2 <= 1024 & floor(tgd2) == tgd2);
            if ~isempty(convInds)
                tgd2(convInds) = constmon.vote.twosComp2dec(dec2bin(tgd2(convInds),10))*1e-9*0.1;
            end

            backwardsScalingInds = find(abs(tgd2)> 1e8);
            if ~isempty(backwardsScalingInds)
                tgd2(backwardsScalingInds) = tgd2(backwardsScalingInds)*(1e-9*0.1)^2;
            end
            
            Eph(34,:) = tgd2;
            Eph(35,:) = tom;
            Eph(36,:) = fit_int;
            
            if suglFlag
                suglInfo1 = data{38};
                suglInfo2 = data{39};
                
                Eph(37,:) = suglInfo1;
                Eph(38,:) = suglInfo2;
            end
            
            if 1; %LSB_RECOVERY
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
            
            
        case 6 % SBAS (RINEX marker: "S")
            
            if (~constellations.SBAS.enabled), continue, end
            
            svprn = data{2};
            year  = data{3};
            year(year<80) = year(year<80)+2000; % See RINEX 3 section 6.10
            year(year<=99) = year(year<=99)+1900; % See RINEX 3 section 6.10
            
            
            month  = data{4};
            day    = data{5};
            hour   = data{6};
            minute = data{7};
            second = data{8};
            
            aGf0 = data{9};
            aGf1 = data{10};
            tom  = data{11}; % See RINEX 3.4, section 6.13
            % If time of message is invalid (garbage data input) use 2 hours
            % before ephemeris time. NOTE: text in 6.13 says to reduce ToM
            % by 604800 (possibly going negative) to refer to same week as
            % ToE refers to... would that be more correct than (toe-7200),
            % as done here?
            if abs(tom) > 86400*7
                tom = toe-7200;
            end
            
            svposX          = data{12};
            svvelXdot       = data{13};
            svaccXdoubledot = data{14};
            svhealth        = data{15};
            
            svposY          = data{16};
            svvelYdot       = data{17};
            svaccYdoubledot = data{18};
            ura             = data{19};
            
            svposZ          = data{20};
            svvelZdot       = data{21};
            svaccZdoubledot = data{22};
            IODN            = data{23};
       
            Eph(1,:)  = svprn;
            Eph(2,:)  = year-2000;
            Eph(3,:)  = month;
            Eph(4,:)  = day;
            Eph(5,:)  = hour;
            Eph(6,:)  = minute;
            Eph(7,:)  = second;
            Eph(8,:)  = aGf0;
            Eph(9,:)  = aGf1;
            Eph(10,:) = tom;
            Eph(11,:) = svposX;
            Eph(12,:) = svvelXdot;
            Eph(13,:) = svaccXdoubledot;
            Eph(14,:) = svhealth;
            Eph(15,:) = svposY;
            Eph(16,:) = svvelYdot;
            Eph(17,:) = svaccYdoubledot;
            Eph(18,:) = ura;
            Eph(19,:) = svposZ;
            Eph(20,:) = svvelZdot;
            Eph(21,:) = svaccZdoubledot;
            Eph(22,:) = IODN;
            
    end % switch constIdx

end % for constIdx

