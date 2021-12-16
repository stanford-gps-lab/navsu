function [obsGnss, dcbCorr] = preprocessGnssObs(obsGnssRaw,corrData,varargin)

p = inputParser;

p.addParameter('epochStart',-Inf);
p.addParameter('epochEnd',Inf);
p.addParameter('downsampleFactor',1);
p.addParameter('obsDesired',[]);
% p.addParameter('duration',Inf);

% parse the results
parse(p, varargin{:});
res = p.Results;
epochStart       = res.epochStart;       % Minimum time of observations
epochEnd         = res.epochEnd;         % Maximum time of observations
downsampleFactor = res.downsampleFactor; % FActor by which to downsample
obsDes           = res.obsDesired;       % RINEX observaiton code cell array... fairly complicated, sorry

%%
obsInds    = repmat([1 2 3 4], 1, 3); % 1 = code, 2 = carrier, 3 = snr, 4 = doppler
signalInds = kron(1:3, ones(1, 4));


%% Extract which measurements should be used
if isempty(obsDes) || 1
    % No specific observations were specified- check whats's available in
    % the input observations

    nConst = 5; % max number of constellations
    nSigs = 3; % max number of signals per constellation
    
    sigList = repmat({'9B'}, nConst, nSigs);

    for cdx = 1:nConst
        switch cdx
            case 1
                % GPS
                sigRankings = {{'1C' '1W' '1Y'};...  % L1
                               {'2W' '2Y' '2P' '2X' '2S' '2L' '5X' '5I' '5Q'}};  % L2/5
            case 2
                % GLONASS
                sigRankings = {{'1C' '1P'};...
                               {'2C' '2P'}};
            case 3
                % GALILEO
                sigRankings = {{'1B' '1C' '1X' '1Z'};...
                               {'5X' '8X' '6X' '7I' '7X' '7Q'}};
            case 4
                % BDS
                sigRankings = {{'1I' '1Q' '1X' '2I' '2Q' '2X' };...
                               {'7I' '6I' '7Q' '6Q' '7X' '6X'}};
                
            otherwise
                continue
        end
        
        indsConst = find(obsGnssRaw.constInds == cdx);
        
        if isempty(indsConst)
            % this constellation isn't even available in the raw data- keep
            % going
            continue;
        end
        
        
        % search for an L1 and up to two other signals to use
        for freqNum = 1:2

            sigNum = freqNum;

            for idx = 1:length(sigRankings{freqNum})
                
                % Get the name of the signal we're looking for
                sigi = sigRankings{freqNum}{idx};

                % check if any measurements available
                if ~isfield(obsGnssRaw.meas, ['C' sigi]) ...
                    || ~any(obsGnssRaw.meas.(['C' sigi])(indsConst, :), 'all')
                    continue;
                end
                
                sigList(cdx, sigNum) = {sigi};
                
                if freqNum == 1
                    break; % only want one signal on L1
                else
                    sigNum = sigNum + 1;
                end
                if sigNum > nSigs
                    break % that's enough signals
                end
            end

        end
    end
    
    % Build the obsDes matrix... which is somewhat wonky
    % Contains Code, Carrier, CN0, Doppler for each signal
    obsDes = repmat({{'C0S'}}, nConst, nSigs*4);
    
    for rdx = 1:nConst
        for idx = 1:nSigs
            obsDes(rdx,((idx-1)*4+1):(idx*4)) = {{['C' sigList{rdx,idx}]} ...
                {['L' sigList{rdx,idx}]} {['S' sigList{rdx,idx}]}  {['D' sigList{rdx,idx}]}};
            
        end
    end
    
    % Just setting a default
%     obsDes  = {{'C1C'} {'L1C'} {'S1C'}  {'D1C'} {'C2S'} {'L2S'} {'S2S'} {'D2S'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
%         {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
%         {'C1B'} {'L1B'} {'S1B'} {'D1B'} {'C7I'} {'L7I'} {'S7I'} {'D7I'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
%         {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
%         {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}  };
    
    
end


%% extract separate measurements
epochs    = obsGnssRaw.epochs;
prns      = obsGnssRaw.PRN(:)'; % make sure it's a row vector
constInds = obsGnssRaw.constInds(:)'; % make sure it's a row vector

nPrn = length(prns);

jdu = unique(floor(navsu.time.epochs2jd(epochs)-0.5)+0.5);
[dayList,YearList] = navsu.time.jd2doy(jdu);

nEpochs = length(epochs);

nObs = size(obsDes,2);
obsOut = zeros(nObs, nEpochs, nPrn);
obsTypes = cell(nObs, nPrn);
for pdx = 1:nPrn
    consti = constInds(pdx);
    if consti > nConst
        % illegal constellation
        warning(['Skipping ', navsu.svprn.convertConstIndName(consti), ...
                 ' satellite ', num2str(prns(pdx)), ...
                 ' (constellation not supported).']);
        continue
    end
    for odx = 1:nObs
        obsSeti = obsDes{consti,odx};
        for osdx = 1:length(obsSeti)
            obsTypei = obsSeti{osdx};
            if ~isfield(obsGnssRaw.meas, obsTypei)
                continue
            end
            if ~any(obsGnssRaw.meas.(obsTypei)(pdx,:))
                continue
            end
            % Save it off.
            obsOut(odx,:,pdx) = obsGnssRaw.meas.(obsTypei)(pdx,:);
            obsTypes{odx,pdx} = obsTypei;
            
            % Don't need to look for additional data
            break
        end
    end
end

% extract frequency number
freqDes2 = nan(size(obsTypes));
haveObs = ~cellfun(@isempty, obsTypes);
freqDes2(haveObs) = cellfun(@(x) str2double(x(2)), obsTypes(haveObs));

[freqs, freqInds] = navsu.svprn.mapSignalFreq(freqDes2(obsInds == 1 | obsInds == 2,:)',prns,constInds,jdu(1));

[freqsDop, ~] = navsu.svprn.mapSignalFreq(freqDes2(obsInds == 4,:)',prns,constInds,jdu(1));


prph12   = squeeze(obsOut(obsInds == 1 | obsInds == 2,:,:));
prphType = squeeze(obsTypes(obsInds == 1 | obsInds == 2,:));
prphSig  = signalInds(obsInds == 1 | obsInds == 2);
prphInd  = obsInds(obsInds == 1 | obsInds == 2);

snr12    = squeeze(obsOut(obsInds == 3 ,:,:));
snrType  = squeeze(obsTypes(obsInds == 3,:));
snrSig   = signalInds(obsInds == 3);

dop12    = squeeze(obsOut(obsInds == 4,:,:));
dopType  = squeeze(obsTypes(obsInds == 4,:));
dopSig   = signalInds(obsInds == 4);


% convert carrier phase from cycles to meters
prph12(prphInd == 2, :, :) = prph12(prphInd == 2, :, :) ...
                             ./ permute(freqs(:, prphInd == 2), [2 3 1]) ...
                             * navsu.constants.c;

% convert doppler from cycles/second to meters/second
dop12 = dop12 ./ permute(freqsDop, [2 3 1]) * navsu.constants.c;


%% Load dcb data
YearListDcb = YearList;
dayListDcb = dayList;
epochDcb = navsu.time.jd2epochs(navsu.time.doy2jd(YearListDcb,dayListDcb))+100;

dcbData = corrData.dcb;
% dcbType = dcbData.type;
% dcbType = 2;


%% Correct L1C-L1P for ISC (using GPS and Galileo MGEX precise products, 
%% this should be the only necessary change)

% Edit Fabian Rothmaier 08/2021: not sure this section is doing anything at
% all, not working with example. dcbType descriptions don't match other
% files.

%     dcbData2 = [];
dcbCorr = zeros(size(prphType));

if ~isempty(dcbData) && dcbData.type == 1 
    % CODE
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if consti == 1
                % GPS
                statCode = [];
            elseif consti == 2
                % GLONASS
                statCode = {'AJAC'};
            else
                % other consts currently not supported
                continue
            end
            
            dcbCorr(idx,jdx) = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if consti == 2 && ~isempty(statCode) && any(strcmp(obsi{1}, {'C1P', 'C2P'}))
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx) + biasi;
            end
            
        end
    end
    
elseif ~isempty(dcbData) && dcbData.type == 3
    % these are relative corrections that need to be further referenced to
    % the L1P-L2P combination (need the TGD term)
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if isempty(obsi{1}) || ~strcmp(obsi{1}(1),'C')
                continue;
            end
            
            if  consti == 1
                % Pull tgd term
            
                % Correct to either C1W or C2W then apply TGD or gamma*TGD
                % from there
                
                bias1 = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1W'},epochDcb,dcbData);
                bias2 = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C2W'},epochDcb,dcbData);
                
                tgd = (bias1-bias2);
                gamma12 = (1575.42/1227.6)^2;
                
                switch obsi{1}
                    case 'C1C'
                        % C1C goes directly to C1W
                        biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1W'},epochDcb,dcbData) ...
                              + tgd;
                    case {'C5Q' 'C5X'}
                        % L5 refers to C1C- the C1C reference must then be
                        % moved to C1W then the tgd can be applied
                        biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1W'},epochDcb,dcbData) ...
                              - navsu.readfiles.findDcbElement(prni,consti,{'C1C'},obsi,epochDcb,dcbData) ...
                              + tgd;
                        
                    case 'C2W'
                        biasi = gamma12*tgd;
                    case 'C2S'
                        biasi = gamma12*tgd ...
                              - navsu.readfiles.findDcbElement(prni,consti,{'C2W'},{'C2S'},epochDcb,dcbData);
                    otherwise
                        biasi = 0;
                end
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if consti == 2
                % All GLONASS DCBs refer to C1C, so just need TGD term for
                % that.
                bias1 = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData);
                bias2 = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C2P'},epochDcb,dcbData);
                
                tgd = (bias1-bias2);
                
                switch obsi{1}
                    case 'C1C'
                        biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData) ...
                              + tgd;
                    case 'C2P'
                        biasi = gamma12*tgd;
                    case 'C2C'
                        biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData) ...
                              - navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C2C'},epochDcb,dcbData) ...
                              + tgd;
                    otherwise
                        biasi = 0;
                end
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            
            if consti == 3
                % All galileo
                'fdaf';
%                 bias1 = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData);
%                 bias2 = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C2P'},epochDcb,dcbData);
%                 
%                 tgd = (bias1-bias2);
                
%                 switch obsi{1}
%                     case 'C1C'
%                         biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData);
%                         biasi = biasi+tgd;
%                     case 'C2P'
%                         biasi = gamma12*tgd;
%                     case 'C2C'
%                         biasi = -navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C2C'},epochDcb,dcbData);
%                         biasi = biasi + navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData);
%                         biasi = biasi+tgd;
%                         
%                     otherwise
%                         biasi = 0;
%                 end
%                 
%                 dcbCorr(idx,jdx) = biasi;
            end
            
        end
    end
    
    
elseif ~isempty(dcbData) && dcbData.type == 5
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if consti == 1 && any(strcmp(obsi{1}, {'C1C', 'C2S'}))
                dcbCorr(idx,jdx) = navsu.readfiles.findDcbElement( ...
                    prni, consti, obsi, {'ABS'}, epochDcb, dcbData);
            end
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if consti == 2 && any(strcmp(obsi{1}, {'C1P', 'C2P'})) && ~isempty(statCode)
                dcbCorr(idx,jdx) = navsu.readfiles.findDcbElement( ...
                    prni, consti, obsi, {'ABS'}, epochDcb, dcbData, statCode);
            end
        end
    end
    
    
elseif ~isempty(dcbData) && dcbData.type == 2
    % STANFORD
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if ~isempty(obsi{1})
                if consti == 1 && strcmp(obsi(1),'C2S')
                    obsi(1) = {'C2X'};
                end
                
                if consti == 2
                    'fdafad';
                end
                                
                dcbCorr(idx,jdx) = navsu.readfiles.findDcbElement(prni,consti,obsi(1),{'ABS'},epochDcb,dcbData);
            end
        end
    end
    
    
elseif ~isempty(dcbData) && dcbData.type == 6
    epochDcb = NaN;
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if isempty(obsi{1})
                continue;
            end
            
            freqi = str2double(obsi{1}(2));
            
            if  consti == 1
                switch freqi
                    case 1
                        if strcmp(obsi{1},'C1W')
                            biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                        else
                            biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData) ...
                                  - navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'ABS'},epochDcb,dcbData);
                        end
                    case 2
                        if strcmp(obsi{1},'C2W')
                            biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                        else
                            biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData) ...
                                  - navsu.readfiles.findDcbElement(prni,consti,{'C2X'},{'ABS'},epochDcb,dcbData);
                        end
                end
                
            elseif consti == 2
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,{'AJAC'});
            end
            
            dcbCorr(idx,jdx) = biasi;
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if any(strcmp(obsi{1}, {'C1P', 'C2P'})) && consti == 2 && ~isempty(statCode)
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx) + biasi;
            end
            
        end
    end
    
    
else
    disp('No DCB''s available :(')
end

    
dcbCorr(isnan(dcbCorr)) = 0;
% Apply the corrections to the observations
prph12i = prph12 - permute(dcbCorr*navsu.constants.c, [1 3 2]);
prph12i(prph12 == 0) = 0;
prph12i(isnan(prph12i)) = 0;


%% Add dual frequency measurements

% signal pairs to include as iono-free combinations
ifPairs = [1 3;
           1 2];

for idx = 1:size(ifPairs,1)
    freq1 = freqs(:,prphSig == ifPairs(idx,1) & prphInd == 1)';
    freq2 = freqs(:,prphSig == ifPairs(idx,2) & prphInd == 1)';
    freqInd1 = freqInds(:,prphSig == ifPairs(idx,1) & prphInd == 1);
    freqInd2 = freqInds(:,prphSig == ifPairs(idx,2) & prphInd == 1);
    
    pr1 = squeeze(prph12i(prphSig == ifPairs(idx,1) & prphInd == 1,:,:));
    pr2 = squeeze(prph12i(prphSig == ifPairs(idx,2) & prphInd == 1,:,:));
    pr1(pr1 == 0) = nan;
    pr2(pr2 == 0) = nan;
    ph1 = squeeze(prph12i(prphSig == ifPairs(idx,1) & prphInd == 2,:,:));
    ph2 = squeeze(prph12i(prphSig == ifPairs(idx,2) & prphInd == 2,:,:));
    ph1(ph1 == 0) = nan;
    ph2(ph2 == 0) = nan;
    
    prifi = (pr1.*freq1.^2-pr2.*freq2.^2)./(freq1.^2-freq2.^2);
    phifi = (ph1.*freq1.^2-ph2.*freq2.^2)./(freq1.^2-freq2.^2);
    freqif = (freq1.^2-freq2.^2)./(freq1-freq2);
    
    prphSig = [prphSig 100*ifPairs(idx,1)+ifPairs(idx,2) 100*ifPairs(idx,1)+ifPairs(idx,2)];
    prphInd = [prphInd 1 2];
    
    prph12i = cat(1,prph12i,permute(prifi, [3 1 2]));
    prph12i = cat(1,prph12i,permute(phifi, [3 1 2]));
    
    freqs = [freqs freqif' freqif'];
    freqInds = [freqInds freqInd1*10+freqInd2 freqInd1*10+freqInd2];
    
    typeAdd = strcat(prphType(prphSig == ifPairs(idx,1) & prphInd == 1,:), ...
                     prphType(prphSig == ifPairs(idx,2) & prphInd == 1,:));
    typeAdd(cellfun(@length,typeAdd) <= 3) = {[]};
    prphType = [prphType; typeAdd];
    
    typeAdd = strcat(prphType(prphSig == ifPairs(idx,1) & prphInd == 2,:), ...
                     prphType(prphSig == ifPairs(idx,2) & prphInd == 2,:));
    typeAdd(cellfun(@length,typeAdd) <= 3) = {[]};
    prphType = [prphType; typeAdd];
end

prph12i(isnan(prph12i)) = 0;
dop12(isnan(dop12)) = 0;
snr12(isnan(snr12)) = 0;

%% Arrange the lock time information (cycle slip detection) for easier use
if ~isempty(obsGnssRaw.tLock)
    lockTime = nan(size(prph12i));
    % Loop through each satellite
    for pdx = 1:length(prns)
        indSat = obsGnssRaw.PRN == prns(pdx) ...
               & obsGnssRaw.constInds == constInds(pdx);
        % Loop through each signal
        for jdx = 1:size(prph12i,1)
            %             indSignal
            sigNamei = prphType{jdx,pdx};
            
            if isempty(sigNamei) %| ~strcmp(sigNamei(1),'L')
                % no obs available- skip this
                continue;
            end
            
            if length(sigNamei) == 3
                % single frequency
                lockTimei = obsGnssRaw.tLock.(sigNamei)(indSat, :);
            else
                % dual frequency- take the minimum of both
                lockTimei = min([obsGnssRaw.tLock.(sigNamei(1:3))(indSat, :);
                                 obsGnssRaw.tLock.(sigNamei(4:6))(indSat, :)]);
            end
            
            lockTime(jdx,:,pdx) = lockTimei;
        end
    end
else
    lockTime = [];
end



%% Put everything into the output structure
obsGnss.PRN       = prns;
obsGnss.constInds = constInds;
obsGnss.epochs    = epochs;
obsGnss.freqs     = freqs';

obsGnss.range.obs         = permute(prph12i, [1 3 2]);
obsGnss.range.rnxCode     = prphType;
obsGnss.range.ind         = repmat(prphInd', 1, size(obsGnss.range.obs,2));
obsGnss.range.sig         = repmat(prphSig', 1, size(obsGnss.range.obs,2));
obsGnss.range.freqs       = freqs';
obsGnss.range.PRN         = repmat(prns,      size(obsGnss.range.obs,1), 1);
obsGnss.range.constInds   = repmat(constInds, size(obsGnss.range.obs,1), 1);
obsGnss.range.lockTime    = permute(lockTime,[1 3 2]);


temp = repelem(navsu.internal.MeasEnum.Code, size(obsGnss.range.ind,1), size(obsGnss.range.ind,2));
temp(obsGnss.range.ind == 1) = navsu.internal.MeasEnum.Code;
temp(obsGnss.range.ind == 2) = navsu.internal.MeasEnum.Carrier;
obsGnss.range.subtype     = temp;

obsGnss.range.ID = navsu.internal.MeasIdGnss(obsGnss.range.PRN, ...
                                             obsGnss.range.constInds, ...
                                             obsGnss.range.sig, ...
                                             temp);

obsGnss.doppler.obs       = permute(dop12,[1 3 2]);
obsGnss.doppler.rnxCode   = dopType;
obsGnss.doppler.sig       = repmat(dopSig', 1, size(obsGnss.doppler.obs,2));
obsGnss.doppler.freqs     = freqsDop';
obsGnss.doppler.PRN       = repmat(prns,      size(obsGnss.doppler.obs,1), 1);
obsGnss.doppler.constInds = repmat(constInds, size(obsGnss.doppler.obs,1), 1);
temp = repelem(navsu.internal.MeasEnum.Doppler, size(obsGnss.doppler.sig,1), size(obsGnss.doppler.sig,2));
obsGnss.doppler.subtype     = temp;

obsGnss.doppler.ID = navsu.internal.MeasIdGnss(obsGnss.doppler.PRN, ...
                                               obsGnss.doppler.constInds, ...
                                               obsGnss.doppler.sig, ...
                                               temp);

obsGnss.snr.obs           = permute(snr12,[1 3 2]);
obsGnss.snr.rnxCode       = repmat(snrType', 1, size(obsGnss.snr.obs,2));
obsGnss.snr.sig           = repmat(snrSig', 1, size(obsGnss.snr.obs,2));
obsGnss.snr.PRN           = repmat(prns, size(obsGnss.snr.obs,1), 1);
obsGnss.snr.constInds     = repmat(constInds, size(obsGnss.snr.obs,1), 1);

%% Pull off only the data we need
%% Pull off data we need
indsGnssMeas = find(obsGnss.epochs >= epochStart & obsGnss.epochs < epochEnd);
indsGnssMeas = indsGnssMeas(1:downsampleFactor:end);
obsGnss.epochs      = obsGnss.epochs(indsGnssMeas);
obsGnss.range.obs   = obsGnss.range.obs(:,:,indsGnssMeas);
obsGnss.doppler.obs = obsGnss.doppler.obs(:,:,indsGnssMeas);
obsGnss.snr.obs     = obsGnss.snr.obs(:,:,indsGnssMeas);

obsGnss.type = navsu.internal.MeasEnum.GNSS;
