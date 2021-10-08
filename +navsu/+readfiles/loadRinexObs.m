function [obsStruc, constellations, time, date, pos, interval, antoff, antmod,...
    rxmod] = loadRinexObs(filename, varargin)
% loadRinexObs
% DESCRIPTION:
%   Parses RINEX observation files.
% INPUT:
%   filename = RINEX observation file
% OPTIONAL INPUTS:
%   'constellations' = constellation object created by
%                    navsu.readfiles.initConstellation and can be used to
%                    select what constellations are parsed
%   'forceRnx3Codes' = boolean- true indicates that for even RINEX version 2
%                    inputs, RINEX 3 observations codes (i.e. C1C) should be
%                    used.
%   'headerOnly'     = boolean - true indicates that only the header should
%                    be read and returned.
%
% OUTPUT:
%   obsStruc = structure containing the parsed observations.  Each field
%            the structure is a different RINEX observation code with the
%            size nPrn x nEpochs
%   constellations = constellation structure containing indexing
%            information about what constellations were parsed and where
%            they are indexed within obsStruc
%   time     = GPS epoch = GPSweek*604800+TOW of each measurement- length
%            nEpochs
%   pos      = rover approximate position pulled from header
%   interval = observation time interval [s]
%   antoff   = antenna offset [m]
%   antmod   = antenna model [string]
%   rxmod    = receiver model
%
% See also:  navsu.readfiles.loadRinexNav
%
% This has been heavily modified from the original goGPS code.
%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
% Portions of code contributed by Damiano Triglione and Stefano Caldera
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



if iscell(filename)
    % If the input was a list of filenames, loop through this code recursively
    % and add each file
    for idx = 1:length(filename)
        [obsStruci, constellationsi, timei, datei, posi, intervali, antoffi, antmodi,...
            rxmodi] = navsu.readfiles.loadRinexObs(filename{idx}, varargin{:});
        
        if idx == 1
            % Initialize outputs
            obsStruc = obsStruci;
            constellations = constellationsi;
            time = timei;
            date = datei;
            pos = posi;
            interval = intervali;
            antoff = antoffi;
            antmod = antmodi;
            rxmod = rxmodi;
        else
            % Concatenate to list of outputs
            
            % Build a new constellations that represents any constellations
            % that have been included at all
            constellationFull = navsu.readfiles.initConstellation(...
                (constellations.GPS.enabled | constellationsi.GPS.enabled),...
                (constellations.GLONASS.enabled | constellationsi.GLONASS.enabled),...
                (constellations.Galileo.enabled | constellationsi.Galileo.enabled),...
                (constellations.BeiDou.enabled | constellationsi.BeiDou.enabled),...
                (constellations.QZSS.enabled | constellationsi.QZSS.enabled),...
                (constellations.SBAS.enabled | constellationsi.SBAS.enabled));
            
            % Time indices for
            nSvFull = length(constellationFull.PRN);
            nEpochsFull = length(time)+length(timei);
            
            % Base observation matrix
            obsEmpty0 = zeros(nSvFull,nEpochsFull);
            
            % Map time indices from individual matrices to the combined one
            tinds1 = 1:length(time);
            tinds2 = (length(time)+1):nEpochsFull;
            
            [~,sinds1] = ismember([constellations.PRN' constellations.constInds'],...
                [constellationFull.PRN' constellationFull.constInds'],'rows');
            [~,sinds2] = ismember([constellationsi.PRN' constellationsi.constInds'],...
                [constellationFull.PRN' constellationFull.constInds'],'rows');
            
            % Loop through each observation type and add it
            obsTypes = unique([fields(obsStruc); fields(obsStruci)]);
            
            obsStrucTemp = [];
            for odx = 1:length(obsTypes)
                obsStrucTemp.(obsTypes{odx}) = obsEmpty0;
                
                if isfield(obsStruc,obsTypes{odx})
                    % Add the old data
                    obsStrucTemp.(obsTypes{odx})(sinds1,tinds1) = ...
                        obsStruc.(obsTypes{odx});
                end
                
                if isfield(obsStruci,obsTypes{odx})
                    % Add the new data
                     % Add the old data
                    obsStrucTemp.(obsTypes{odx})(sinds2,tinds2) = ...
                        obsStruci.(obsTypes{odx});
                end
            end
            % The temporary structure is now the current one
            obsStruc = obsStrucTemp;
            
            % Collect the rest of the data
            constellations = constellationFull;
            time = [time; timei];
            date = [date; datei];
            pos = [pos posi];
            interval = [interval intervali];
            antoff = [antoff antoffi];
            antmod = [antmod antmodi];
            rxmod = [rxmod rxmodi];
        end
    end
    
else
   
    % Parse a single file
    p = inputParser;
    p.addParameter('headerOnly', false);
    p.addParameter('forceRnx3Codes', false);
    p.addParameter('constellations',[]);
    
    % parse the results
    parse(p, varargin{:});
    res = p.Results;
    headerOnly = res.headerOnly;
    forceRnx3Codes = res.forceRnx3Codes;
    constellations = res.constellations;
    
    %variable initialization
    nEpochsAdd = 3000;
    nEpochs = 3000;
    time = NaN(nEpochs,1);
    tow = NaN(nEpochs,1);
    week = NaN(nEpochs,1);
    date = NaN(nEpochs,6);
    pos = zeros(3,1);
    interval = zeros(1,1);
    antoff = zeros(3,1);
    antmod = cell(1,1);
    rxmod = cell(1,1);
    
    %open RINEX observation file
    fid = fopen(filename,'r');
    
    %parse RINEX header
    [obs_type, pos(:,1), basic_info, interval(1,1), sysId, antoff(:,1), ...
        antmod{1,1}, ~, rxmod{1,1}, ~] = navsu.readfiles.rinexParseHeader(fid);
    
    %check the availability of basic data to parse the RINEX file
    if (basic_info == 0)
        error(['RINEX file ' filename ': basic data is missing in the file header'])
    end
    
    % Pull each signal name out
    [obsColumns, nObsTypes, obsTypes] = navsu.readfiles.rinexFindObsType(obs_type, sysId);
    
    % Add the LLI field for each carrier-phase frequency if this is RINEX 2
    if any(find(contains(obsTypes,'L1'))) &&  isempty(sysId)
        obsTypes{end+1} = 'LLI1';
        obsColumns{end+1} = 'LLI1';
    end
    if any(find(contains(obsTypes,'L2')))  && isempty(sysId)
        obsTypes{end+1} = 'LLI2';
        obsColumns{end+1} = 'LLI2';
    end
    if any(find(contains(obsTypes,'L5'))) && isempty(sysId)
        obsTypes{end+1} = 'LLI5';
        obsColumns{end+1} = 'LLI5';
    end
    
    % Convert to RINEX 3 Codes if necessary and desired
    if isempty(sysId) && forceRnx3Codes
        constTypes = {'G','R','E','C','J','S'};
        
        obsTypes3 = {};
        for cdx = 1:length(constTypes)
            obsColumnsi = navsu.readfiles.convertRinex3ObsCodes(obsColumns,constTypes{cdx});
            obsColumns3.(constTypes{cdx}) = obsColumnsi;
            
            obsTypes3 = unique([obsTypes3; obsColumnsi(~cellfun(@isempty,obsColumnsi))]);
        end
        
        obsColumns = obsColumns3;
        obsTypes = obsTypes3;
        sysId = 'CONVERTED_3';
    end
    
    if isempty(sysId) % RINEX v2.xx
        obsColumnsMat = repmat(1:length(obsTypes),6,1);
    else % RINEX v3.xx
        constTypes = {'G','R','E','C','J','S'};
        obsColumnsMat = zeros(length(constTypes),length(obsTypes));
        for cdx = 1:length(constTypes)
            for odx = 1:length(obsTypes)
                if isfield(obsColumns,constTypes{cdx})
                    %                 coli = obsColumns.(constTypes{cdx}).(obsTypes{odx});
                    coli = find(~cellfun(@isempty,strfind(obsColumns.(constTypes{cdx}) , obsTypes{odx})));
                    
                    if ~isempty(coli)
                        obsColumnsMat(cdx,odx) = coli(1);
                    end
                end
            end
        end
    end
    
    % If specific constellation usage hasn't been specified, set it based on
    % what's available.
    if isempty(constellations)   && ~isempty(sysId) && ~any(strcmp(sysId,'CONVERTED_3'))
        % RINEX 3
        constellations = navsu.readfiles.initConstellation(...
            ismember({'G'},sysId),  ismember({'R'},sysId) , ...
            ismember({'E'},sysId),  ismember({'C'},sysId), ...
            ismember({'J'},sysId),  ismember({'S'},sysId));
        
    elseif isempty(constellations) && (isempty(sysId) || any(strcmp(sysId,'CONVERTED_3')))
        warning('Defaulting to GPS and GLONASS only for RINEX 2.  Consider using a ''constellation'' input')
        % RINEX 2 - just default to using GPS GLONASS only
        constellations = navsu.readfiles.initConstellation(1,1,0,0,0,0);
    end
    
    nSatTot = constellations.nEnabledSat;
    
    % initialize storage variables
    for odx = 1:length(obsTypes)
        obsStruc.(obsTypes{odx}) = NaN(nSatTot,nEpochs);
    end
    obsOut = nan(nSatTot, nEpochs, length(obsTypes));
    
    if headerOnly
        % if we only want the header, ignore the rest of the file
        fclose(fid);
    else
        % Parse the rest of the file
        k = 1;
        while (~feof(fid))
            %read data for the current epoch (ROVER)
            [time(k,1), date(k,:), num_sat, sat, sat_types, tow(k,1)] = navsu.readfiles.rinexGetEpoch(fid);
            
            if (k > nEpochs)
                obsOutTemp = NaN(size(obsOut, 1), ...
                                 nEpochs + nEpochsAdd, ...
                                 size(obsOut, 3), ...
                                 size(obsOut, 4));
                obsOutTemp(:, 1:size(obsOut, 2), :, :) = obsOut;
                obsOut = obsOutTemp;
                
                dateTemp = NaN(nEpochs + nEpochsAdd, ...
                               size(date, 2), ...
                               size(date, 3));
                dateTemp(1:size(date, 1), :, :) = date;
                date = dateTemp;
                
                towTemp = NaN(nEpochs+nEpochsAdd,size(tow,2),size(tow,3));
                towTemp(1:size(tow, 1), :) = tow;
                tow = towTemp;
                
                timeTemp = NaN(nEpochs+nEpochsAdd,size(time,2),size(time,3));
                timeTemp(1:size(time, 1), :) = time;
                time = timeTemp;
                
                weekTemp = zeros(nEpochs+nEpochsAdd,size(week,2),size(week,3));
                weekTemp(1:size(week, 1), :) = week;
                week = weekTemp;
                
                nEpochs = nEpochs  + nEpochsAdd;
            end
            
            %read ROVER observations
            [~,obsMati] = navsu.readfiles.rinexGetObs(fid, num_sat, sat, sat_types, obsColumns, ...
                nObsTypes, constellations, obsColumnsMat, obsTypes);
            
            obsOut(:,k,:) = obsMati;
            
            k = k + 1;
        end
        
        %GPS week number
        week(:,1) = navsu.time.epochs2gps(navsu.time.cal2epochs(date(:,:)));
        
        %observation rate
        if (interval(:,1) == 0)
            interval(:,1) = round((median(time(2:k-1,1) - time(1:k-2,1)))*1000)/1000;
        end
        
        %close RINEX file
        fclose(fid);
    end
    
    for odx = 1:length(obsTypes)
        obsStruc.(obsTypes{odx}) = squeeze(obsOut(:,:,odx));
    end
    
    [time, date, obsStruc, interval] = ...
        navsu.readfiles.rinexSyncObs(time, week, date, obsStruc, interval);
end

for odx = 1:length(obsTypes)
    obsStruc.(obsTypes{odx}) = squeeze(obsOut(:,:,odx));
end

% Peel off LLI flags from .0001's digit (relies on modified rinexGetObs observation masking hack!)
obsFields = fieldnames(obsStruc);
idxObsTypesCarrier = ... % LLI applies to carrier observations only (obs codes 'Lxx')
    find(strcmp(cellfun(@(x) x(1), obsFields, 'UniformOutput', false),'L'));
for idxLLI = 1:length(idxObsTypesCarrier)
    dummy = obsStruc.(obsFields{idxObsTypesCarrier(idxLLI)});
    idxBlankFlag = (dummy>0 & dummy<10); % handle blank lines (which parse as zeros) with non-zero LLI flags
    dummy(idxBlankFlag) = dummy(idxBlankFlag)/10000; % shift flag into .0001's digit position
    flags = round(10*rem(dummy*1000,1),1);
    obsTypes{end+1} = ['LLI' obsFields{idxObsTypesCarrier(idxLLI)}]; %#ok<AGROW>
    obsStruc.(obsTypes{end}) = flags;
    obsStruc.(obsFields{idxObsTypesCarrier(idxLLI)}) = round(dummy, 3); % strip off flag digit from observations themselves to restore original values
end

[time, date, obsStruc, interval] = ...
    navsu.readfiles.rinexSyncObs(time, week, date, obsStruc, interval);

end

