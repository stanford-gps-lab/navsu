function [dcbData, tecData] = readIonex(filename)
% readIonex
% DESCRIPTION:
%   File parser for IGS IONEX files. 
% INPUT:
%   filename    - the name of the file (along with the path) to parse.
% OUTPUT:
%   dcbData     - structure containing DCB data- similar to those from 
%                 0navsu.readfiles.parseDcbBsxFile
%   tecData     - structure containing ionospheric data
%
% See also: navsu.readfiles.loadDcb, navsu.svOrbitClock

%
% CURRENTLY ONLY PARSING SATELLITE DCB DATA.. SORRY!

fid = fopen(filename,'r');
dcbData = [];
if fid == -1
    return;
end

while ~feof(fid)
    % Just loop through until beginning of section of interest
    linetxt = fgetl(fid);
    
    if length(linetxt) < 14
        continue;
    end
    
    if contains(linetxt,'EPOCH OF FIRST MAP')
        datai = textscan(linetxt,'%f %f %f %f %f %f');
        
        startJdi = navsu.time.cal2jd(datai{1},datai{2},datai{3}+datai{4}/24+datai{5}/(24*60)+datai{6}/(24*60*60));
        [startDyi,startYri] = navsu.time.jd2doy(startJdi);
        startEpochi = navsu.time.jd2epochs(startJdi);
    end
    
    if contains(linetxt,'EPOCH OF LAST MAP')
        datai = textscan(linetxt,'%f %f %f %f %f %f');
        
        endJdi = navsu.time.cal2jd(datai{1},datai{2},datai{3}+datai{4}/24+datai{5}/(24*60)+datai{6}/(24*60*60));
        [endDyi,endYri] = navsu.time.jd2doy(endJdi);
        endEpochi = navsu.time.jd2epochs(endJdi);
    end
    
    if contains(linetxt,'HGT1 / HGT2 / DHGT')
        datai = textscan(linetxt,'%f %f %f');
        
        hgt1 = datai{1};
        hgt2 = datai{2};
        dhgt = datai{3};
    end
    if contains(linetxt,'LAT1 / LAT2 / DLAT')
        datai = textscan(linetxt,'%f %f %f');
        lat1 = datai{1};
        lat2 = datai{2};
        dlat = datai{3};
    end
    if contains(linetxt,'LON1 / LON2 / DLON')
        datai = textscan(linetxt,'%f %f %f');
        
        lon1 = datai{1};
        lon2 = datai{2};
        dlon = datai{3};
    end
    
    if contains(linetxt,'# OF MAPS IN FILE')
        datai = textscan(linetxt,'%f');
        
        nMaps = datai{1};
    end
    
    % Offset between geodetic marker and antenna reference point
    if contains(linetxt, 'START OF AUX DATA') ...
        && contains(upper(linetxt), 'DIFFERENTIAL CODE BIASES')
        % Check if we have a signal definition
        if strcmp(linetxt(27:31), '     ')
            obs1i = 'C1W';
            obs2i = 'C2W';
        else
%             pair = upper(linetxt(27:31));
%             pair(strfind(pair,'P')) = 'W';
            
            obs1i = ['C' linetxt(28:-1:27)];
            obs2i = ['C' linetxt(31:-1:31)];
        end
               
        nAllocInit = 1e2;
        nAlloc = nAllocInit;
        %         pairData = nan(nAllocInit,4);
        SVNs  = zeros(nAllocInit,1);
        PRNs  = zeros(nAllocInit,1);
        sites = cell(nAllocInit,1);
        domes = cell(nAllocInit,1);
        obs1  = cell(nAllocInit,1);
        obs2  = cell(nAllocInit,1);
        startYr    = zeros(nAllocInit,1);
        startDy    = zeros(nAllocInit,1);
        startEpoch = zeros(nAllocInit,1);
        endYr      = zeros(nAllocInit,1);
        endDy      = zeros(nAllocInit,1);
        endEpoch   = zeros(nAllocInit,1);
        bias       = zeros(nAllocInit,1);
        stdDev     = zeros(nAllocInit,1);
        satFlag    = zeros(nAllocInit,1);
        const      = cell(nAllocInit,1);
        constInd   = zeros(nAllocInit,1);
        
        ind = 1;
        % get lines until end of section
        while ~contains(linetxt, 'END OF AUX DATA')
            
            % get next line
            linetxt = fgetl(fid);
            
            if contains(upper(linetxt), 'COMMENT') ...
               && strcmpi(linetxt(1:34), 'Reference observables for GPS    :')
                obs1i = linetxt(36:38);
                obs2i = linetxt(40:42);
            end
            
            if contains(upper(linetxt),'PRN / BIAS / RMS')
                datai = textscan(linetxt,'%1s %2f %f %f',1);
                
                if ind > nAlloc
                    SVNs  = [SVNs; zeros(nAllocInit,1)];
                    PRNs  = [PRNs; zeros(nAllocInit,1)];
                    sites = [sites; cell(nAllocInit,1)];
                    domes = [domes; cell(nAllocInit,1)];
                    obs1  = [obs1; cell(nAllocInit,1)];
                    obs2  = [obs2; cell(nAllocInit,1)];
                    startYr    = [startYr; zeros(nAllocInit,1)];
                    startDy    = [startDy; zeros(nAllocInit,1)];
                    startEpoch = [startEpoch; zeros(nAllocInit,1)];
                    endYr      = [endYr; zeros(nAllocInit,1)];
                    endDy      = [endDy; zeros(nAllocInit,1)];
                    endEpoch   = [endEpoch; zeros(nAllocInit,1)];
                    bias       = [bias; zeros(nAllocInit,1)];
                    stdDev     = [stdDev; zeros(nAllocInit,1)];
                    satFlag    = [satFlag; zeros(nAllocInit,1)];
                    const      = [const; cell(nAllocInit,1)];
                    constInd   = [constInd; zeros(nAllocInit,1)];
                    nAlloc = nAlloc+nAllocInit;
                end
                
                consti  = datai{1}{1};
                prni    = datai{2};
                biasi   = datai{3};
                stdDevi = datai{4};
                
                % get constellation index
                constIndi = navsu.svprn.convertConstIndName(consti,1);
                
                SVNs(ind) = navsu.svprn.prn2svn(prni,startJdi,constIndi);
                PRNs(ind) = prni;
                sites{ind} = NaN;
                domes{ind} = NaN;
                obs1{ind} = obs1i;
                obs2{ind} = obs2i;
                startYr(ind) = startYri;
                startDy(ind) = startDyi;
                startEpoch(ind) = startEpochi;
                endYr(ind) = endYri;
                endDy(ind) = endDyi;
                endEpoch(ind) = endEpochi;
                bias(ind) = biasi;
                stdDev(ind) = stdDevi;
                satFlag(ind) = 1;
                const{ind} = consti;
                constInd(ind) = constIndi;
                
                ind = ind+1;
            end
        end
        
        % Save everything in a struct!
        dcbData.SVNs        = SVNs(1:ind-1);
        dcbData.PRNs        = PRNs(1:ind-1);
        dcbData.sites       = sites(1:ind-1);
        dcbData.domes       = domes(1:ind-1);
        dcbData.obs1        = obs1(1:ind-1);
        dcbData.obs2        = obs2(1:ind-1);
        dcbData.startYr     = startYr(1:ind-1);
        dcbData.startDy     = startDy(1:ind-1);
        dcbData.startEpoch  = startEpoch(1:ind-1);
        dcbData.endYr       = endYr(1:ind-1);
        dcbData.endDy       = endDy(1:ind-1);
        dcbData.endEpoch    = endEpoch(1:ind-1);
        dcbData.bias        = bias(1:ind-1);
        dcbData.stdDev      = stdDev(1:ind-1);
        dcbData.satFlag     = satFlag(1:ind-1);
        dcbData.const       = const(1:ind-1);
        dcbData.constInd    = constInd(1:ind-1);
    end
    
    % once we reach the end of the header, initialize stuff for the TEC
    % map.
    if contains(linetxt,'END OF HEADER')
        % initialize stuff
        
        latVec = lat1:dlat:lat2;
        lonVec = lon1:dlon:lon2;
        if dhgt == 0
            hgtVec= hgt1;
        else
            hgtVec = hgt1:dhgt:hgt2;
        end
        
        tecMap = nan(length(latVec),length(lonVec),length(hgtVec),nMaps);
        
        epochs = nan(nMaps,1);
        datevecs = nan(nMaps,6);
       
        
    end
    
    % also parse the actual iono data
    if contains(linetxt,'START OF TEC MAP')
        
        datai = textscan(linetxt,'%f');
        mapi = datai{1};
        
        % Get epoch of current map
        linetxt = fgetl(fid);
        datai = textscan(linetxt,'%f %f %f %f %f %f');
        datevecs(mapi,:) = cell2mat(datai);
        epochs(mapi) = navsu.time.cal2epochs(datevecs(mapi,1),datevecs(mapi,2),...
            datevecs(mapi,3),datevecs(mapi,4),datevecs(mapi,5),datevecs(mapi,6));
        
        tecMapi = nan(length(latVec),length(lonVec),length(hgtVec));
        
        linetxt = fgetl(fid);
        firstLine = 1;
        while ~contains(linetxt,'END OF TEC MAP')
            if contains(linetxt,'LAT/LON1/LON2/DLON/H')
                if ~firstLine
                   % add this entry to the map 
                   latInd   = latVec == lati;
                   hgtInd   = hgtVec == hgti;
                   
                   tecMapi(latInd, :, hgtInd) = mapEntryi;
                end
                
                % Initialize 
                datai = cell2mat(textscan(linetxt,'%f %f %f%f %f'));
                lati = datai(1);
%                 lonsi = datai(2):datai(4):datai(3);
                hgti   = datai(5);
                
                mapEntryi = [];
                firstLine = 0;
            
            else 
                datai = textscan(linetxt,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ');
                
                mapEntryi = [mapEntryi [datai{:}]];
                
            end
            linetxt = fgetl(fid);
        end
        % Need to do thte last entry
        latInd   = latVec == lati;
        hgtInd   = hgtVec == hgti;
        tecMapi(latInd, :, hgtInd) = mapEntryi;
        
        % Put it all in the big map!
        tecMap(:,:,:,mapi) = tecMapi;
    end
    
end

tecData.tecMap = tecMap;
tecData.latVec = latVec;
tecData.lonVec = lonVec;
tecData.epochs = epochs;
tecData.datevecs = datevecs;
tecData.hgtVec  = hgtVec;

fclose(fid);

end














