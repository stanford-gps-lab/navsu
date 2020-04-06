function pvt = ReadSP3Mixed(filename, cmToApcFlag,strictConstNumFlag,constellationOut,atxData)
% TrueEphem = readNGA( filename )
% TrueEphem is a n x 12 matrix. Each row of TrueEphem includes
% [ GPS week, GPS sec, PRN, x, y, z, Clock drift, dx, dy, dz, d(Clock drift), Event Flag ]
% No mathematical processing --- the value of each element is identical to that in the NGA file

% Written by: Liang Heng  05/28/2009

% functions called: err_chk, utc2leap

pvt = [];
% if (nargin <= 1)
%     outputFormat = [];
% end
[fid, message] = fopen(filename, 'rt');
if (fid == -1)
    fprintf(2, 'Error open %s: %s\n', filename, message);
    return
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

NumSV = -1;
% a = [];
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
            %         case '++'
            %             if length(a) < NumSV
            %                 a = [a; sscanf(tline(10:end), '%d')];
            %             end
        case '* '
            break
    end
end
% a = a(1 : NumSV);
consts = 'GRECJ';
C = [];
while ~feof(fid)
    tline = fgetl(fid);
    if ~strcmp(tline(1),'*') && ~strcmp(tline,'EOF')
        const = tline(2);
        
        constNum = strfind(consts,const);
        
        %%%% CHANGING THIS 3/6/2018 FOR JAX DATA
%         Ci = textscan(tline(3:end),'%n%n%n%n%n');
%         dataTemp = [Ci{:}];

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
        %         NumSV = NumSV;
    end
end
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
    'Event',array(:, 5) .* 0); % SHOULD ACTUALLY CHECK FOR EVENTS

pvt.Event(abs(pvt.clock_bias) >= 0.999999) = true;
pvt.Event(any(pvt.position == 0, 2)) = true;

%% Correct APC
if cmToApcFlag
    % Pull sun position for each time in SP3 file
    
    if ~contains(path,'\mice\src\mic')
        navsu.thirdparty.attachToMice();
    end
    t = GPS_seconds:Epoch_interval:(GPS_seconds+NumEpochs*Epoch_interval);
    jd = gps2jd(GPS_week_num,t);
    sunpos = zeros(3,length(t));
    for i = 1:length(jd)
        et       = cspice_str2et(['jd ' num2str(jd(1))]);
        sunposi     = cspice_spkezr( 'sun',et , 'itrf93', 'none', 'earth');
        sunpos(:,i) = sunposi(1:3)*1000;
    end
    
    PRNs = unique(array(:,1));
    
    if ~isempty(atxData)
        constNum = constellationOut;
        atxDatai= atxData([atxData.type] == constNum);
        epochi = GPS_seconds+604800*GPS_week_num;
    end
    
    for pdx = 1:NumSV
        prn = PRNs(pdx);
        if ~isempty(atxData)
            svni = prn2svn(prn,jd(1),constNum);
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
        else
            if constellationOut == 1
                offset = antennaOffsetGPS(svn2block(prn2svn(prn,jd)));
                %             offset = [0 0 0]';
            elseif constellationOut == 2
                offset = apcOffsetGlonass(prn,GPS_seconds+GPS_week_num*604800);
            elseif constellationOut == 3
                offset = apcOffsetGalileo(prn);
            elseif constellationOut == 4
                offset = [0.6 0 1.1]';
            else
                offset = [0 0 0]';
            end
        end
        
        inds = find(array(:,1) == prn);
        
        if any(any(~isnan(pvt.position(inds,:))))
            for tdx = 1:NumEpochs
                sunposi = sunpos(:,tdx);
                svPosi  = pvt.position(inds(tdx),:)';
                
                % Build body axis rotation matrix
                e = (sunposi-svPosi)./norm(sunposi-svPosi);
                k = -svPosi./norm(svPosi);
                % yhat = cross(e,k)./norm(cross(e,k));
                j = cross(k,e)/norm(cross(k,e));
                i = cross(j,k)/norm(cross(j,k));
                
                R = [i j k];
                offsetECEF = R*offset;
                pvt.position(inds(tdx),:) = svPosi+offsetECEF;
                
            end
        end
    end
end



























