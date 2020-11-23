function [solEstimate,statEcc,solEpochs,siteRcvr,siteAtx] = parseSscFile(filename)

% Parser for IGS station solution files


fid = fopen(filename,'r');

% stationNames= {};
% stationPositions = [];
solEstimate = [];
statEcc     = [];
solEpochs   = [];
siteRcvr    = [];
siteAtx     = [];


while ~feof(fid)
    % Just loop through until beginning of section of interest
    linetxt = fgetl(fid);
    
%     if strcmp(filename(end-2:end),'ssc')
%         indsComp = 1:length(linetext);
%     else
%        indsComp = 1:34; 
%     end
%     
    % Offset between geodetic marker and antenna reference point
    if strcmp('+SITE/ECCENTRICITY',linetxt(1:(min([length('+SITE/ECCENTRICITY') length(linetxt)]))))
        % Start reading!
        % Header line
        linetxt = fgetl(fid);
        
        stationNamesEcc = {};
        stationEcc = [];
        % first line
        linetxt = fgetl(fid);
        % get lines until end of section
        while ~strcmp('-SITE/ECCENTRICITY',linetxt(1:(min([length('-SITE/ECCENTRICITY') length(linetxt)]))))
            % Read data from current line
            stationNamei = linetxt(2:5);
            
            stationNamesEcc = [stationNamesEcc; {stationNamei}];
            
            % Save epochs and reference and offset
            startYr = str2double(linetxt(17:18));
            % convert to actual year
            if startYr > 80
                startYr = startYr + 1900;
            else
                startYr = startYr + 2000;
            end
            startDy = str2double(linetxt(20:22));
            startSc = str2double(linetxt(24:28));
            
            endYr  = str2double(linetxt(30:31));
            if endYr > 80
                endYr = endYr + 1900;
            else
                endYr = endYr + 2000;
            end
            endDy  = str2double(linetxt(33:35));
            endSc  = str2double(linetxt(37:41));
            
            ref = linetxt(43:45);
            refFlag = strcmp(ref,'UNE');
            
            dxU = str2double(linetxt(47:54));
            dxN = str2double(linetxt(56:63));
            dxE = str2double(linetxt(65:72));
            
            stationEcc = [stationEcc; [startYr startDy startSc endYr endDy endSc refFlag dxU dxN dxE]];
            
            % get next line
            linetxt = fgetl(fid);
            
        end
        statEcc.name = stationNamesEcc;
        statEcc.data = stationEcc;
        
    end
    
    if strcmp('+SOLUTION/EPOCHS',linetxt(1:(min([length('+SOLUTION/EPOCHS') length(linetxt)]))))
        % Start reading!
        % Header line
        linetxt = fgetl(fid);
        
        solEpochData = [];
        solEpochNames = [];
        nSolEpochs = 1;
        
        % first line
        linetxt = fgetl(fid);
        % get lines until end of section
        while ~strcmp('-SOLUTION/EPOCHS',linetxt(1:(min([length('-SOLUTION/EPOCHS') length(linetxt)])))) 
            % Read data from current line
            stationNamei = linetxt(2:5);
            solEpochNames = [solEpochNames; {stationNamei}];

            PT = linetxt(7:8);
            solni = str2num(linetxt(10:13));
            
            meanYr = str2num(linetxt(43:44));
            if meanYr > 80
                meanYr = meanYr + 1900;
            else
                meanYr = meanYr + 2000;
            end
            meanDy = str2num(linetxt(46:48));
            meanSc = str2num(linetxt(50:54));
            
            % get next line
            linetxt = fgetl(fid);
            
            solEpochData = [solEpochData; [meanYr meanDy meanSc solni]];
        end
        solEpochs.name = solEpochNames;
        solEpochs.data = solEpochData;
        
    end
    
    
    if strcmp('+SITE/RECEIVER',linetxt(1:(min([length('+SITE/RECEIVER') length(linetxt)])))) 
        % first line
        linetxt = fgetl(fid);
        statCodes     = [];
        startTimes    = [];
        endTimes      = [];
        receiverTypes = [];
        serialNumbers = [];
        firmwares     = [];
        
        if strcmp(linetxt(1),'*')
            % this is just a header line! keep reading.
            linetxt = fgetl(fid);
        end
        
        % get lines until end of section
        while ~strcmp('-SITE/RECEIVER',linetxt(1:(min([length('-SITE/RECEIVER') length(linetxt)]))))
            % Read data from current line
            stationNamei = linetxt(2:5);

            startYear = str2num(linetxt(17:18));
            if startYear > 80
                startYear = startYear + 1900;
            else
                startYear = startYear + 2000;
            end
            startDy = str2num(linetxt(20:22));
            startSc = str2num(linetxt(24:28));
           
            endYear = str2num(linetxt(30:31));
            endDy = str2num(linetxt(33:35));
            endSc = str2num(linetxt(37:41));
            if endYear+endDy+endSc == 0
                endYear = Inf;
            end
            
            if endYear > 80
                endYear = endYear + 1900;
            else
                endYear = endYear + 2000;
            end
            
            
            receiverTypei = linetxt(43:63);
            
            serialNumberi = linetxt(64:69);
            
            firmwarei   = linetxt(70:length(linetxt));
            
            % get next line
            linetxt = fgetl(fid);
            
            % Save data;
            statCodes     = [statCodes; {stationNamei}];
            startTimes    = [startTimes; [startYear startDy startSc]];
            endTimes      = [endTimes; [endYear endDy endSc]];
            receiverTypes = [receiverTypes; {receiverTypei}];
            serialNumbers = [serialNumbers; {serialNumberi}];
            firmwares     = [firmwares; {firmwarei}];
        end
        
        siteRcvr.statCodes     = statCodes;
        siteRcvr.startTimes    = startTimes;
        siteRcvr.endTimes      = endTimes;
        siteRcvr.receiverTypes = receiverTypes;
        siteRcvr.serialNumbers = serialNumbers;
        siteRcvr.firmwares     = firmwares;
    end
    
    if strcmp('+SOLUTION/ESTIMATE',linetxt(1:(min([length('+SOLUTION/ESTIMATE') length(linetxt)])))) 
        % Start reading!
        % Header line
        linetxt = fgetl(fid);
        
        % first line
        linetxt = fgetl(fid);
        stationNames = [];
        
        % get lines until end of section
        while ~strcmp('-SOLUTION/ESTIMATE',linetxt(1:(min([length('-SOLUTION/ESTIMATE') length(linetxt)])))) ...
                && length(linetxt) >= 80
            % Read data from current line
            stationNamei = linetxt(15:18);
            
            component = linetxt(8:11);
            value = str2num(linetxt(48:68));
            
            % Check if station is already on list
            statDx = find(strcmp(stationNames,stationNamei));
            
            if isempty(statDx)
                stationNames = [stationNames; {stationNamei}];
                statDx = length(stationNames);
            end
            
            soln = str2double(linetxt(23:26));
            
            refYr = str2num(linetxt(28:29));
            if refYr > 80
                refYr = refYr + 1900;
            else
                refYr = refYr + 2000;
            end
            refDy = str2num(linetxt(31:33));
            refSc = str2num(linetxt(35:39));
            refEpoch = navsu.time.jd2epochs(navsu.time.doy2jd(refYr,refDy))+refSc;
            
            
            componentList = {'STAX','STAY','STAZ','VELX','VELY','VELZ'};
            compDx = find(strcmp(componentList,component));
            
            
            
            if ~isempty(compDx)
                stationPositions(statDx,compDx) = value;
                stationPositions(statDx,7) = soln;
                stationPositions(statDx,8) = refEpoch;
            end
            
            % get next line
            linetxt = fgetl(fid);
            
        end
        
        solEstimate.name = stationNames;
        solEstimate.data = stationPositions;
%         solEstimate.
        
    end
    
    if strcmp('+SITE/ANTENNA',linetxt(1:(min([length('+SITE/ANTENNA') length(linetxt)]))))
        % first line
        linetxt = fgetl(fid);
        statCodes     = [];
        startTimes    = [];
        endTimes      = [];
        antennaTypes = [];
        serialNumbers = [];
        
        if strcmp(linetxt(1),'*')
            % this is just a header line! keep reading.
            linetxt = fgetl(fid);
        end
        
        % get lines until end of section
        while ~strcmp('-SITE/ANTENNA',linetxt(1:(min([length('-SITE/ANTENNA') length(linetxt)])))) 
            % Read data from current line
            stationNamei = linetxt(2:5);

            startYear = str2num(linetxt(17:18));
            if startYear > 80
                startYear = startYear + 1900;
            else
                startYear = startYear + 2000;
            end
            startDy = str2num(linetxt(20:22));
            startSc = str2num(linetxt(24:28));
           
            endYear = str2num(linetxt(30:31));
            endDy = str2num(linetxt(33:35));
            endSc = str2num(linetxt(37:41));
            if endYear+endDy+endSc == 0
                endYear = Inf;
            end
            
            if endYear > 80
                endYear = endYear + 1900;
            else
                endYear = endYear + 2000;
            end
            
            antennaTypei = linetxt(43:62);
            
            serialNumberi = linetxt(64:end);
            
            % get next line
            linetxt = fgetl(fid);
            
            % Save data;
            statCodes     = [statCodes; {stationNamei}];
            startTimes    = [startTimes; [startYear startDy startSc]];
            endTimes      = [endTimes; [endYear endDy endSc]];
            antennaTypes = [antennaTypes; {antennaTypei}];
            serialNumbers = [serialNumbers; {serialNumberi}];
        end
        
        siteAtx.statCodes     = statCodes;
        siteAtx.startTimes    = startTimes;
        siteAtx.endTimes      = endTimes;
        siteAtx.receiverTypes = antennaTypes;
        siteAtx.serialNumbers = serialNumbers;
    end 
    
end

fclose(fid);

end
