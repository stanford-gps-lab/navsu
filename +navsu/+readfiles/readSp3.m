function pvt = readSp3(filename, cmToApcFlag,strictConstNumFlag,constellationOut,atxData)
% readSp3
% DESCRIPTION:
%   Parses .sp3 precise orbit and clock files!
% INPUT:
%   filename    - name of the .sp3 file to parse
% OPTIONAL INPUTS:
%   cmToApcFlag - default is false.  Flag to indicate to offset the
%                 positions from the default center of mass to antenna
%                 phase center.
%   strictConstNumFlag - default is false.  Flag to indicate to force each
%                 constellation to have a specific default number of satellites.
%                 nGPS = 32, nGlonass = 24, nGalileo = 36, nBeidou = 35
%   constellationOut - default is 1.  Constellation index indicating which
%                 constellation to parse. 12345 = GRECS
%   atxData     - antenna phase center data for converting from center of
%                 mass to antenna phase center.  This should be from an IGS
%                 .atx file
% OUTPUT:
%   pvt           - structure containing all of the parsed data
%    .filename    - name of the file that was parsed
%    .GPS_week_num - vector of GPS weeks corresponding to unique output
%                   epochs
%    .GPS_seconds - vector of GPS times of week corresponding to unique
%                   output epochs
%    .Epoch_interval - nominal spacing between epochs [s]
%    .NumSV       - number of satellites included in output
%    .NumEpochs   - number of unique epochs included in output
%    .PRN         - vector list of PRNs in output
%    .clock_bias  - vector of clock biases in seconds- each element
%                   corresponds to values of .GPS_seconds and .PRN
%    .position    - [length(.PRN) x 3] matrix of satellite positions in m
%    .Event       - vector of detected estimator events
%
% See also:  navsu.svOrbitClock, navsu.readfiles.loadPEph

pvt = [];

[fid, message] = fopen(filename, 'rt');
if (fid == -1)
    fprintf(2, 'Error open %s: %s\n', filename, message);
    return
end
if nargin < 2
    cmToApcFlag = false;
end

if (nargin < 3)
    strictConstNumFlag = false;
end

%GRECJ

if nargin < 4
    constellationOut = 1;
end

if nargin < 5
    atxData = [];
end

if cmToApcFlag && isempty(atxData)
    error(['Antenna phase center data (IGS .atx file) required to ' ...
        'translate from center of mass to antenna phase center']);
end

NumSV = -1;
% Default constellation is GPS
const = 'GPS';
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    switch tline(1:2)
        case '#a'
            PRNformat = ' %d';
            DateTime = sscanf(tline(4:end), '%f' ,6);
            NumEpochs = str2double(tline(33:39));
            if tline(3) == 'P'
                dataFormat = 'P %n%n%n%n%n%s%*[^\n]\n';
            elseif tline(3) == 'V'
                dataFormat = 'P %n%n%n%n%n%s%*[^\n]\nV %*[^\n]\n';
            else
                fprintf(2, 'Error: Unsupported file format %s.\n', tline(1:3));
                fclose(fid);
                return
            end
        case {'#c' '#d'}
            PRNformat = 'G%d';
            DateTime = sscanf(tline(4:end), '%f' ,6);
            NumEpochs = str2double(tline(33:39));
            if tline(3) == 'P'
                dataFormat = 'PG%n%n%n%n%n%n%n%n%n%s%*[^\n]\n';
            elseif tline(3) == 'V'
                dataFormat = 'PG%n%n%n%n%n%n%n%n%n%s%*[^\n]\nVG%*[^\n]\n';
            else
                fprintf(2, 'Error: Unsupported file format %s.\n', tline(1:3));
                fclose(fid);
                return
            end
        case '%c'
            if ~strcmp(tline(4),'c')
                PRNformat = [tline(4) '%d'];
                dataFormat(2) = tline(4);
                if strcmp(tline(4),'R')
                    const = 'GLO';
                end
                if strcmp(tline(4),'M')
                    const = 'MIX';
                end
                
            end
        case '##'
            GPS_week_num = str2double(tline(4:7));
            GPS_seconds = str2double(tline(9:23));
            Epoch_interval = str2double(tline(25:38));
        case '+ '
            if NumSV < 0
                NumSV = str2double(tline(2:8));
                if NumSV < 1
                    fprintf(2, 'Error: The number of SVs is %d.\n', NumSV);
                    fclose(fid);
                    return
                end
                s = sscanf(tline(10:end), PRNformat);
                continue
            end
            if length(s) < NumSV
                s = [s; sscanf(tline(10:end), PRNformat)];
            end
        case '* '
            break
    end
end
consts = 'GRECJ';
C = [];
while ~feof(fid)
    tline = fgetl(fid);
    if ~strcmp(tline(1),'*') && ~strcmp(tline,'EOF')
        const = tline(2);
        
        constNum = strfind(consts,const);
        
        Ci = textscan(tline(3:end),'%n%n%n%n%n%n%n%n%n%n%n%n%n');
        dataTemp = [Ci{1:5}];
        
        if ~isempty(dataTemp)
            C = [C; [constNum dataTemp]];
        end
    end
end

C = C(C(:,1) == constellationOut,2:end);
array = C;
fclose(fid);

if strictConstNumFlag
    % Limit number of SVs to maximum for that constellation (sometimes GLONASS
    % actually has 25-26 SVs) and ensure that each block of time has that
    % number of SVs\
    if constellationOut == 1 % GPS
        outArray = nan(32*NumEpochs,9);
        setInds = [0; find(diff(array(:,1)) < 1); length(array(:,1))];
        for tdx = 1:NumEpochs
            for prn = 1:32
                setIndsi = (setInds(tdx)+1):setInds(tdx+1);
                arrayi = array(setIndsi(array(setIndsi,1) == prn),:);
                if ~isempty(arrayi)
                    outArray((tdx-1)*32+prn,1:length(arrayi)) = arrayi;
                else
                    outArray((tdx-1)*32+prn,:) = [prn nan(1,8)];
                end
            end
        end
        array = outArray;
        NumSV = 32;
    end
    
    if constellationOut == 2 % GLONASS
        outArray = nan(24*NumEpochs,9);
        setInds = [0; find(diff(array(:,1)) < 1); length(array(:,1))];
        for tdx = 1:NumEpochs
            for prn = 1:24
                setIndsi = (setInds(tdx)+1):setInds(tdx+1);
                arrayi = array(setIndsi(array(setIndsi,1) == prn),:);
                if ~isempty(arrayi)
                    outArray((tdx-1)*24+prn,1:length(arrayi)) = arrayi;
                else
                    outArray((tdx-1)*24+prn,:) = [prn nan(1,8)];
                end
            end
        end
        array = outArray;
        NumSV = 24;
    end
    if constellationOut == 3 % Galileo
        outArray = nan(36*NumEpochs,5);
        setInds = [0; find(diff(array(:,1)) < 1); length(array(:,1))];
        for tdx = 1:NumEpochs
            for prn = 1:36
                setIndsi = (setInds(tdx)+1):setInds(tdx+1);
                arrayi = array(setIndsi(array(setIndsi,1) == prn),:);
                if ~isempty(arrayi)
                    outArray((tdx-1)*36+prn,:) = arrayi;
                else
                    outArray((tdx-1)*36+prn,:) = [prn nan(1,4)];
                end
            end
        end
        array = outArray;
        NumSV = 36;
    end
    
    if constellationOut == 4 % Beidou
        NumSV = 35;
        outArray = nan(NumSV*NumEpochs,5);
        setInds = [0; find(diff(array(:,1)) < 1); length(array(:,1))];
        if length(setInds) > 2
            for tdx = 1:NumEpochs
                for prn = 1:NumSV
                    setIndsi = (setInds(tdx)+1):setInds(tdx+1);
                    arrayi = array(setIndsi(array(setIndsi,1) == prn),:);
                    if ~isempty(arrayi)
                        outArray((tdx-1)*NumSV+prn,:) = arrayi;
                    else
                        outArray((tdx-1)*NumSV+prn,:) = [prn nan(1,4)];
                    end
                end
            end
        end
        array = outArray;
    end
end

epochs = ones(NumSV, 1) * (GPS_seconds + ...
    GPS_week_num*7*24*3600 + ...
    Epoch_interval .* (0:NumEpochs-1));
epochs = epochs(:);
constInds = constellationOut*ones(size(epochs));

pvt = struct('filename', filename, ...
    'DateTime', DateTime', ...
    'GPS_week_num', GPS_week_num, ...
    'GPS_seconds', GPS_seconds, ...
    'Epoch_interval', Epoch_interval, ...
    'NumSV', NumSV, ...
    'NumEpochs', NumEpochs, ...
    'PRN', array(:, 1), ...
    'clock_bias', array(:, 5) .* 1e-6, ...
    'position', array(:, 2:4) .* 1000, ...
    'velocity',nan(size(array(:,2:4))), ...
    'Event',array(:, 5) .* 0,...
    'epochs', epochs,...
    'constInds',constInds);

pvt.Event(abs(pvt.clock_bias) >= 0.999999) = true;
pvt.Event(any(pvt.position == 0, 2)) = true;

%% Correct APC
if cmToApcFlag
    % Pull sun position for each time in SP3 file
    
    %     if ~contains(path,'\mice\src\mic')
    %         navsu.thirdparty.attachToMice();
    %     end
    t = GPS_seconds:Epoch_interval:(GPS_seconds+NumEpochs*Epoch_interval);
    jd = navsu.time.gps2jd(GPS_week_num,t);
    epochsSun = navsu.time.gps2epochs(GPS_week_num*ones(size(t)),t);
    %     sunpos = zeros(3,length(t));
    %     for i = 1:length(jd)
    %         et          = cspice_str2et(['jd ' num2str(jd(1))]);
    %         sunposi     = cspice_spkezr( 'sun',et , 'itrf93', 'none', 'earth');
    %         sunpos(:,i) = sunposi(1:3)*1000;
    %     end
    % ECEF sun position
    sunpos = navsu.geo.sunVecEcef(jd)';
    
    
    
    
    PRNs = unique(array(:,1));
    
    if ~isempty(atxData)
        constNum = constellationOut;
        atxDatai= atxData([atxData.type] == constNum);
        epochi = GPS_seconds+604800*GPS_week_num;
    end
    
    for pdx = 1:NumSV
        prn = PRNs(pdx);
        % Use the IGS antenna phase center file
        svni = navsu.svprn.prn2svn(prn,jd(1),constNum);
        if isnan(svni)
            offset = nan(3,1);
        else
            adx = find([atxDatai.svn] == svni & [atxDatai.epochStart] <= epochi & [atxDatai.epochEnd] > epochi);
            
            if ~isempty(adx)
                offset = (atxDatai(adx).apc(1,:)')*1e-3;
            else
                offset = nan(3,1);
            end
        end
        
        inds = find(array(:,1) == prn);
        
        if any(any(~isnan(pvt.position(inds,:))))
            for tdx = 1:NumEpochs
                sunposi = sunpos(:,tdx);
                svPosi  = pvt.position(inds(tdx),:)';
                
                R = navsu.geo.svLocalFrame(svPosi',epochsSun(tdx),sunposi);
                
                offsetECEF = R*offset;
                
                pvt.position(inds(tdx),:) = svPosi+offsetECEF;
            end
        end
    end
end



























