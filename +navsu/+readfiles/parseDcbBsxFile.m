function dcbData = parseDcbBsxFile(filename,skipRx)
% parseDcbBsxFile
% DESCRIPTION:
%   Parses a specified IGS differential code bias .bsx file.  
% INPUT:
%   filename    - name of the file to be parsed (including path)
%   skipRx      - boolean, where TRUE tells the function not to parse all
%                 of the receiver-specific entries (saves time!)
% OUTPUT:
%   dcbData     - structure containing the DCB data
%
% See also: navsu.readfiles.loadDcb, navsu.ftp.download

fid = fopen(filename,'r');
dcbData = [];
if fid == -1
   return; 
end

if nargin < 2
    skipRx = 0;
end

% Default bias mode is relative
biasMode = 'REL';

while ~feof(fid) 
    % Just loop through until beginning of section of interest
    linetxt = fgetl(fid);
        
    if length(linetxt) < 14
        continue;
    end
    
    if strcmp('+BIAS/DESCRIPT',linetxt(1:14))
       linetxt = fgetl(fid);
       
       while ~strcmp('-BIAS/DESCRIPTION',linetxt(1:17))
           % Check for bias mode
           
           if contains(linetxt,'BIAS_MODE')
               if contains(linetxt,'ABSOLUTE')
                  biasMode = 'ABS'; 
               end
           end
           linetxt = fgetl(fid);
       end
    end
    
    % Offset between geodetic marker and antenna reference point
    if strcmp('+BIAS/SOLUTION',linetxt(1:14))
        % Start reading!
        % Header line
        linetxt = fgetl(fid);
        
        areDomes = contains(linetxt,'DOMES____');
        
        if strcmp(linetxt(1),'*')
            % this is just another header line- get the next line to begin
            % parsing
            % first line
            linetxt = fgetl(fid);
        end
        
        % Initialize storage
        nAllocInit = 1e6;
        nAlloc = nAllocInit;
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
        
        if strcmp(biasMode,'REL')
            obs2Format = '%3s';
        else
            obs2Format = [];
        end
        
        
        % get lines until end of section
        while ~strcmp('-BIAS/SOLUTION',linetxt(1:14))
            % Read data from current line
            
            % Check constellation
            consti = linetxt(7);
            constIndi = strfind('GRECSJ',consti);
            
            if contains(linetxt,'IFB') || strcmp(linetxt(1),'*')
                linetxt = fgetl(fid);
                
                continue;
            end
            
            if strcmp(linetxt(16:19),'    ')
                % this is a satellite.
%                 temp = textscan(linetxt,'%3s %1s %3f %1s %2f %3s %3s %2f %1s %3f %1s %5f %2f %1s %3f %1s %5f %2s %f %f');
                temp = textscan(linetxt,['%3s %1s %3f %1s %2f %3s ' obs2Format ' %f %1s %f %1s %f %f %1s %f %1s %f %2s %f %f']);
                
                svni = temp{3};
                satFlagi = 1;
                prni   = temp{5};
                sitesi = '';
                domesi = NaN;
                rxIndOffset = 0;
            elseif strcmp(linetxt(21:29),'         ') || ~areDomes && ~skipRx
               
                satFlagi = 0;
                domesi = NaN;
                
                if strcmp(linetxt(8:10),'   ')
                    temp = textscan(linetxt,['%3s %1s %1s %4s %3s ' obs2Format ' %f %1s %f %1s %f %f %1s %f %1s %f %2s %f %f']);
                    
                    svni    = NaN;
                    prni    = NaN;
                    
                    sitesi = temp{4}{1};
                    rxIndOffset = -1;
                else
                    temp = textscan(linetxt,['%3s %1s %f %1s %f %4s %3s ' obs2Format ' %f %1s %f %1s %f %f %1s %f %1s %f %2s %f %f']);
                    
                    svni    = temp{3};
                    prni    = temp{5};
                    
                    sitesi = temp{6}{1};
                    rxIndOffset = 1;
                    
                end
                
            elseif ~skipRx
                temp = textscan(linetxt,['%3s %1s %1s %4s %9s %3s ' obs2Format '  %2f %1s %3f %1s %5f %2f %1s %3f %1s %5f %2s %f %f']);

                satFlagi = 0;
                svni    = NaN;
                prni    = NaN;
                
                sitesi = temp{4}{1};
                domesi = temp{5}{1};
                
                rxIndOffset = 0;
            else
               linetxt = fgetl(fid);
                continue;
            end
            if ~strcmp(biasMode,'REL')
                % Stick in an element for obs2
                
                temp = [temp(1:6+rxIndOffset) {{'ABS'}} temp(7+rxIndOffset:end)];
                
            end
            
            
            obs1i = temp{6+rxIndOffset}{1};
            obs2i = temp{7+rxIndOffset}{1};
            
%             Save epochs and reference and offset
            startYri = temp{8+rxIndOffset};
            % convert to actual year
            if startYri > 80 && startYri < 1000
                startYri = startYri + 1900;
            elseif startYri < 1000
                startYri = startYri + 2000;
            end
            startDyi = temp{10+rxIndOffset};
            startEpochi = navsu.time.jd2epochs(navsu.time.doy2jd(startYri,startDyi));
            
            %             Save epochs and reference and offset
            endYri = temp{13+rxIndOffset};
            % convert to actual year
            if endYri > 80 && endYri < 1000
                endYri = endYri + 1900;
            elseif endYri < 1000
                endYri = endYri + 2000;
            end
            endDyi = temp{15+rxIndOffset};
            endEpochi = navsu.time.jd2epochs(navsu.time.doy2jd(endYri,endDyi));
            
            
            % Bias
            biasi = temp{19+rxIndOffset};

            % Standard Deviation
            stdDevi = temp{20+rxIndOffset};

            % Save everything
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
            SVNs(ind) = svni;
            PRNs(ind) = prni;
            sites{ind} = sitesi;
            domes{ind} = domesi;
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
            satFlag(ind) = satFlagi;
            const{ind} = consti;
            constInd(ind) = constIndi;
            
            ind = ind+1;
            % get next line
            linetxt = fgetl(fid);
            
        end
        
        % Cut off the extra data at the end that was preallocated
        
        
        
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
    
end

dcbData.biasMode = biasMode;

fclose(fid);

end