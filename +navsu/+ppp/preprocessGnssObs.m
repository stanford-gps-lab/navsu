function [obsGnss, dcbCorr] = preprocessGnssObs(obsGnssRaw,corrData,varargin)

% this is mostly a wrapper for svPosFromProd
p = inputParser;

p.addParameter('epochStart',-Inf);
p.addParameter('epochEnd',Inf);
p.addParameter('downsampleFactor',1);
p.addParameter('obsDesired',[]);

% parse the results
parse(p, varargin{:});
res = p.Results;
epochStart       = res.epochStart;       % Minimum time of observations
epochEnd         = res.epochEnd;         % Maximum time of observations
downsampleFactor = res.downsampleFactor; % FActor by which to downsample
obsDes           = res.obsDesired;       % RINEX observaiton code cell array... fairly complicated, sorry

%%
obsInds    = repmat([1 2 3 4],1,3); % 1 = code, 2 = carrier, 3 = snr, 4 = doppler
signalInds = kron(1:3,ones(1,4));


%%
if isempty(obsDes) || 1
    % No specific observations were specified- check whats's available in
    % the input observations

    sigList = repmat({'9B'},5,3);

    for cdx = 1:5
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
            otherwise
                continue
        end
        
        indsConst = find(obsGnssRaw.constInds == cdx);
        
        if isempty(indsConst)
            % this constellation isn't even available in the raw data- keep
            % going
            continue;
        end
        
        % Search for an L1 signal to use 
        sig1Rankings = sigRankings{1};
        for idx = 1:size(sig1Rankings)
           % Just check for the signal at all
           sigi = sig1Rankings{idx};
           
           if ~isfield(obsGnssRaw.meas,['C' sigi])
               continue;
           end
           % Check if this constellation is actually available for this
           % signal
           
           sigAvail = any(any(obsGnssRaw.meas.(['C' sigi])(indsConst,:)));
           
           if sigAvail
               sigList(cdx,1) = {sigi};
               break;
           end
        end

        % Search for dual frequency signals to use
        sigNum = 2;
        sig1Rankings = sigRankings{2};
        for idx = 1:length(sig1Rankings)
           % Just check for the signal at all
           sigi = sig1Rankings{idx};
           
           if ~isfield(obsGnssRaw.meas,['C' sigi])
               continue;
           end
           % Check if this constellation is actually available for this
           % signal
           
           sigAvail = any(any(obsGnssRaw.meas.(['C' sigi])(indsConst,:)));
           
           if sigAvail
               sigList(cdx,sigNum) = {sigi};
               sigNum = sigNum+1;
           end
           if sigNum == 4
               break;
           end
        end
    end
    
    % Build the obsDes matrix... which is somewhat wonky
    obsDes = repmat({{'C0S'}},5,12);
    
    for rdx = 1:size(sigList,1)
        for idx = 1:size(sigList,2)
            obsDes(rdx,((idx-1)*4+1):(idx*4)) = {{['C' sigList{rdx,idx}]} ...
                {['L' sigList{rdx,idx}]} {['S' sigList{rdx,idx}]}  {['D' sigList{rdx,idx}]}};
            
        end
    end
    
    % Just setting a default
    %           1                               2                               3                               4                                5
%     obsDes  = {{'C1C'} {'L1C'} {'S1C'}  {'D1C'} {'C2S'} {'L2S'} {'S2S'} {'D2S'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
%         {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
%         {'C1B'} {'L1B'} {'S1B'} {'D1B'} {'C7I'} {'L7I'} {'S7I'} {'D7I'} {'C2W'} {'L2W'} {'S2W'} {'D2W'}
%         {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}
%         {'C1C'} {'L1C'} {'S1C'} {'D1C'} {'C2C'} {'L2C'} {'S2C'} {'D2C'} {'C2P'} {'L2P'} {'S2P'} {'D2P'}  };
    
    
end


%%
epochs    = obsGnssRaw.epochs;
prns      = obsGnssRaw.PRN;
constInds = obsGnssRaw.constInds;

nPrn = length(prns);

jdu = unique(floor(navsu.time.epochs2jd(epochs)-0.5)+0.5);
[dayList,YearList] = navsu.time.jd2doy(jdu);

nEpochs = length(epochs);

nObs = size(obsDes,2);
obsOut = zeros(nObs,nEpochs,nPrn);
obsTypes = cell(nObs,nPrn);
for pdx = 1:nPrn
    prni = prns(pdx);
    consti = constInds(pdx);
    for odx = 1:nObs
        obsSeti = obsDes{consti,odx};
        for osdx = 1:length(obsSeti)
            obsTypei = obsSeti{osdx};
            if ~isfield(obsGnssRaw.meas,obsTypei)
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

freqDes2 = nan(size(obsTypes));
for idx = 1:size(freqDes2,1)
    for jdx = 1:size(freqDes2,2)
        if ~isempty(obsTypes{idx,jdx})
            freqDes2(idx,jdx) = str2double(obsTypes{idx,jdx}(2));
        end
    end
end
[freqs, freqInds] = navsu.svprn.mapSignalFreq(freqDes2(obsInds == 1 | obsInds == 2,:)',prns,constInds,jdu(1));

[freqsDop, freqIndsDop] = navsu.svprn.mapSignalFreq(freqDes2(obsInds == 4,:)',prns,constInds,jdu(1));


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

c = navsu.constants.c;

% convert carrier phase from cycles to meters
for idx = 1:size(prph12,1)
    if prphInd(idx) == 2
        prph12(idx,:,:) = squeeze(prph12(idx,:,:)).*c./freqs(:,idx)';
    end
end
% convert doppler from cycles/second to meters/second
for idx = 1:size(dop12,1)
    dop12(idx,:,:) = squeeze(dop12(idx,:,:)).*c./freqs(:,min(find(prphSig == dopSig(idx))))';
end


%% Load dcb data
YearListDcb = YearList;
dayListDcb = dayList;
epochDcb = navsu.time.jd2epochs(navsu.time.doy2jd(YearListDcb,dayListDcb))+100;

dcbData = corrData.dcb;
dcbType = corrData.settings.dcbSource;

%% Correct L1C-L1P for ISC (using GPS and Galileo MGEX precise products, this should be the only necessary change)
%     dcbData2 = [];
dcbCorr = zeros(size(prphType));

if dcbType == 4 && ~isempty(corrData.dcb)
    % CODE
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if  consti == 1
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if consti == 2
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,{'AJAC'});
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if strcmp(obsi{1},'C1P') && consti == 2 && ~isempty(statCode)
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
            
            if strcmp(obsi{1},'C2P') && consti == 2 && ~isempty(statCode)
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
        end
    end
    
    dcbCorr = settings.dcbUse*dcbCorr*c;
    
    dcbCorr(isnan(dcbCorr)) = 0;
    % Apply the corrections to the observations
    prph12i = prph12-permute(repmat(dcbCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
    
elseif dcbType == 2 && ~isempty(corrData.dcb)
    % DLR
    % these are relative corrections that need to be further referenced to
    % the L1P-L2P combination (need the TGD term)
    dcbCorr = zeros(size(prphType));
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if isempty(obsi{1}) | ~strcmp(obsi{1}(1),'C')
                continue;
            end
            freqi = str2num(obsi{1}(2));
            
            if  consti == 1
                switch obsi{1}
                    case 'C1C'
                        biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1W'},epochDcb,dcbData);
                    case 'C2W'
                        biasi = 0;
                    case 'C2S'
                        biasi = -navsu.readfiles.findDcbElement(prni,consti,{'C2W'},{'C2S'},epochDcb,dcbData);
                    otherwise
                        biasi = 0;
                end
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if consti == 2
                switch obsi{1}
                    case 'C1C'
                        biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C1P'},epochDcb,dcbData);
                    case 'C2P'
                        biasi = 0;
                    case 'C2C'
                        biasi = -navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'C2C'},epochDcb,dcbData);
                    otherwise
                        biasi = 0;
                end
                
                dcbCorr(idx,jdx) = biasi;
            end
            
        end
    end
    
    
    dcbCorr = dcbCorr*c;
    
    dcbCorr(isnan(dcbCorr)) = 0;
    
    % Apply the corrections to the observations
    prph12i = prph12-permute(repmat(dcbCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
    
elseif dcbType == 5 && ~isempty(corrData.dcb)
    dcbCorr = zeros(size(prphType));
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if strcmp(obsi{1},'C1C') && consti == 1
                %                     biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'C1W'},epochDcb,dcbData);
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if strcmp(obsi{1},'C2S') && consti == 1
                %                     biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'C1W'},epochDcb,dcbData);
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if strcmp(obsi{1},'C1P') && consti == 2 && ~isempty(statCode)
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
            
            if strcmp(obsi{1},'C2P') && consti == 2 && ~isempty(statCode)
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
        end
    end
    
    
    dcbCorr(isnan(dcbCorr)) = 0;
    % Apply the corrections to the observations
    prph12i = prph12-c*permute(repmat(dcbUse*dcbCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
    
elseif dcbType == 2 && ~isempty(corrData.dcb)
    % STANFORD
    dcbCorr = zeros(size(prphType));
    
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
                
                biasi = c*navsu.readfiles.findDcbElement(prni,consti,obsi(1),{'ABS'},epochDcb,dcbData);
                
                dcbCorr(idx,jdx) = biasi;
            end
        end
    end
    
    dcbCorr(isnan(dcbCorr)) = 0;
    fullCorr = dcbCorr;
    
    prph12i = prph12-permute(repmat(fullCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
elseif dcbType == 6 && ~isempty(corrData.dcb)
    epochDcb = NaN;
    
    for idx = 1:size(prphType,1)
        for jdx = 1:size(prphType,2)
            prni = prns(jdx);
            consti = constInds(jdx);
            obsi = prphType(idx,jdx);
            
            if isempty(obsi{1})
                continue;
            end
            
            freqi = str2num(obsi{1}(2));
            
            if  consti == 1
                switch freqi
                    case 1
                        if strcmp(obsi{1},'C1W')
                            biasi = navsu.readfiles.findDcbElement(prni,consti,{'C1W'},{'ABS'},epochDcb,dcbData);
                        else
                            biasi = -navsu.readfiles.findDcbElement(prni,consti,{'C1C'},{'ABS'},epochDcb,dcbData)+...
                                navsu.readfiles.findDcbElement(prni,consti,{'C1W'},{'ABS'},epochDcb,dcbData);
                        end
                    case 2
                        if strcmp(obsi{1},'C2W')
                            biasi = navsu.readfiles.findDcbElement(prni,consti,{'C2W'},{'ABS'},epochDcb,dcbData);
                        else
                            biasi = -navsu.readfiles.findDcbElement(prni,consti,{'C2X'},{'ABS'},epochDcb,dcbData)+...
                                navsu.readfiles.findDcbElement(prni,consti,{'C2W'},{'ABS'},epochDcb,dcbData);
                        end
                end
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            if consti == 2
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,{'AJAC'});
                
                dcbCorr(idx,jdx) = biasi;
            end
            
            % GLONASS adjustments for IGS stations (if using CODE IFB/DCB)
            if strcmp(obsi{1},'C1P') && consti == 2 && ~isempty(statCode)
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
            
            if strcmp(obsi{1},'C2P') && consti == 2 && ~isempty(statCode)
                biasi = navsu.readfiles.findDcbElement(prni,consti,obsi,{'ABS'},epochDcb,dcbData,statCode);
                dcbCorr(idx,jdx) = dcbCorr(idx,jdx)+biasi;
            end
        end
    end
    
    dcbCorr = dcbCorr*c;
    
    dcbCorr(isnan(dcbCorr)) = 0;
    % Apply the corrections to the observations
    prph12i = prph12-permute(repmat(dcbCorr,1,1,size(prph12,2)),[1 3 2]);
    prph12i(prph12 == 0) = 0;
    prph12i(isnan(prph12i)) = 0;
    
else
    disp('No DCB''s available :(')
    prph12i = prph12;
end

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
    
    typeAdd = strcat(prphType(prphSig == ifPairs(idx,1) & prphInd == 1,:),prphType(prphSig == ifPairs(idx,2) & prphInd == 1,:));
    typeAdd(cellfun(@length,typeAdd) <= 3) = {[]};
    prphType = [prphType; typeAdd];
    
    typeAdd = strcat(prphType(prphSig == ifPairs(idx,1) & prphInd == 2,:),prphType(prphSig == ifPairs(idx,2) & prphInd == 2,:));
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
        indSat = find(obsGnssRaw.PRN == prns(pdx) & obsGnssRaw.constInds == constInds(pdx));
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
                lockTimei = obsGnssRaw.tLock.(sigNamei)(indSat,:);
            else
                % dual frequency- take the minimum of both
                lockTimei = min([obsGnssRaw.tLock.(sigNamei(1:3))(indSat,:);
                    obsGnssRaw.tLock.(sigNamei(4:6))(indSat,:)]);
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

obsGnss.range.obs         = permute(prph12i,[1 3 2]);
obsGnss.range.rnxCode     = prphType;
obsGnss.range.ind         = repmat(prphInd',1,size(obsGnss.range.obs,2));
obsGnss.range.sig         = repmat(prphSig',1,size(obsGnss.range.obs,2));
obsGnss.range.freqs       = freqs';
obsGnss.range.PRN         = repmat(prns,size(obsGnss.range.obs,1),1);
obsGnss.range.constInds   = repmat(constInds,size(obsGnss.range.obs,1),1);
obsGnss.range.lockTime    = permute(lockTime,[1 3 2]);

obsGnss.doppler.obs       = permute(dop12,[1 3 2]);
obsGnss.doppler.rnxCode   = repmat(dopType',1,size(obsGnss.doppler.obs,2));
obsGnss.doppler.sig       = repmat(dopSig',1,size(obsGnss.doppler.obs,2));
obsGnss.doppler.freqs     = freqsDop';
obsGnss.doppler.PRN       = repmat(prns,size(obsGnss.doppler.obs,1),1);
obsGnss.doppler.constInds = repmat(constInds,size(obsGnss.doppler.obs,1),1);

obsGnss.snr.obs           = permute(snr12,[1 3 2]);
obsGnss.snr.rnxCode       = repmat(snrType',1,size(obsGnss.snr.obs,2));
obsGnss.snr.sig           = repmat(snrSig',1,size(obsGnss.snr.obs,2));
obsGnss.snr.PRN           = repmat(prns,size(obsGnss.snr.obs,1),1);
obsGnss.snr.constInds     = repmat(constInds,size(obsGnss.snr.obs,1),1);

%% Pull off only the data we need
%% Pull off data we need
indsGnssMeas = find(obsGnss.epochs >= epochStart & obsGnss.epochs < epochEnd);
indsGnssMeas = indsGnssMeas(1:downsampleFactor:end);
obsGnss.epochs      = obsGnss.epochs(indsGnssMeas);
obsGnss.range.obs   = obsGnss.range.obs(:,:,indsGnssMeas);
obsGnss.doppler.obs = obsGnss.doppler.obs(:,:,indsGnssMeas);
obsGnss.snr.obs     = obsGnss.snr.obs(:,:,indsGnssMeas);

obsGnss.type = navsu.internal.MeasEnum.GNSS;

%% Now, convert to a cell array output :)
% outCells = cell(length(obsGnss.epochs),1);
% cell0 = obsGnss;
% cell0.epochs = [];
% cell0.range.obs = [];
% cell0.doppler.obs = [];
% cell0.snr.obs = [];
% cell0.type = navsu.internal.MeasEnum.GNSS;
%
% for idx = 1:length(obsGnss.epochs)
%     celli = cell0;
%     celli.epochs = obsGnss.epochs(idx);
%     celli.range.obs = obsGnss.range.obs(:,:,idx);
%     celli.doppler.obs = obsGnss.doppler.obs(:,:,idx);
%     celli.snr.obs     = obsGnss.snr.obs(:,:,idx);
%
%     if ~isempty(obsGnss.range.lockTime)
%        obsGnss.range.lockTime = obsGnss.range.lockTime(:,:,idx);
%     end
%
%     outCells{idx} = celli;
% end

% 'dafa';
