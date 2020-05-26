function outData = readSPANdata(filename,varargin)

% file extension
[~,~,fileExt] = fileparts(filename);

% Parse varargin
parser = inputParser;
parser.addParameter('RAWIMUSXA', true);
parser.addParameter('RAWIMUA',  true);
parser.addParameter('RANGEA',   true);
parser.addParameter('VEHICLEDATASA', true);
parser.addParameter('constellations',[]);
parser.addParameter('correctHalfCycle',false);
parser.parse(varargin{:});
res = parser.Results;

RAWIMUSXA = res.RAWIMUSXA;
RAWIMUA   = res.RAWIMUA;
RANGEA    = res.RANGEA;
VEHICLEDATASA = res.VEHICLEDATASA;
constellations = res.constellations;
correctHalfCycle = res.correctHalfCycle;

textdata = importdata(filename);

% Vehicle data log :)
if VEHICLEDATASA && ~isempty(find(contains(textdata,'%VEHICLEDATASA'), 1))
    indsi  = find(contains(textdata,'%VEHICLEDATASA'));
    
    formati = '%f%f%f%f%f%f%f%f%s%s%s%f%f%f%f%f%f%s';

    tempi = textdata(indsi);
    
    parseOut = cellfun(@(x) textscan(x,formati,'Delimiter',',;'),cellfun(@(x) x(16:end-9),tempi,...
        'Uni',0),'Uni',0);
    
    parseOut = cat(1,parseOut{:});
    
    outData.VEHICLEDATASA.headerWeek    = cell2mat(parseOut(:,1));
    outData.VEHICLEDATASA.headerTow     = round(cell2mat(parseOut(:,2)));
    outData.VEHICLEDATASA.SpeedFrontLeft  = cell2mat(parseOut(:,3));
    outData.VEHICLEDATASA.SpeedFrontRight = cell2mat(parseOut(:,4));
    outData.VEHICLEDATASA.SpeedRearLeft   = cell2mat(parseOut(:,5));
    outData.VEHICLEDATASA.SpeedRearRight  = cell2mat(parseOut(:,6));
    outData.VEHICLEDATASA.SteeringAngle = cell2mat(parseOut(:,7));
    outData.VEHICLEDATASA.VehicleSpeed  = cell2mat(parseOut(:,8));
    outData.VEHICLEDATASA.Stationary         = cat(1,parseOut{:,9});
    outData.VEHICLEDATASA.TransmissionSetting = cat(1,parseOut{:,10});
    outData.VEHICLEDATASA.ParkingBrakeStatus = cat(1,parseOut{:,11});
    
    outData.VEHICLEDATASA.epochs  = navsu.time.gps2epochs(...
        outData.VEHICLEDATASA.headerWeek,outData.VEHICLEDATASA.headerTow);

else
   outData.VEHICLEDATASA = struct('headerWeek',[],'headerTow',[],...
       'SpeedFrontLeft',[],'SpeedFrontRight',[],'SpeedRearLeft',[],...
       'SpeedRearRight',[],'SteeringAngle',[],'VehicleSpeed',[],...
       'Stationary',[],'TransmissionSetting',[],'ParkingBrakeStatus',[],...
       'epochs',[]);
end



if RAWIMUSXA && ~isempty(find(contains(textdata,'%RAWIMUSXA'), 1))
    % find imudata
    indsi  = find(contains(textdata,'%RAWIMUSXA'));
    formati = '%f,%f;%f,%f,%f,%f,%8s,%f,%f,%f,%f,%f,%f';
    
    tempi = textdata(indsi);
    
    parseOut = cellfun(@(x) textscan(x,formati),cellfun(@(x) x(12:end-9),tempi,'Uni',0),'Uni',0);
    
    % The seventh column is the imu status and is a string- pull it out
    % here
    data1 = cellfun(@(x) x(1:6),parseOut,'UniformOutput',false);
    data2 = cellfun(@(x) x(8:end),parseOut,'UniformOutput',false);
    
    imuTable = cell2mat([cat(1,data1{:}) repelem({0},length(data1),1) cat(1,data2{:}) ]);
    
    imuStatus = cellfun(@(x) x(7),parseOut,'UniformOutput',false);
    
    outData.RAWIMUSXA.headerWeek = imuTable(:,1);
    outData.RAWIMUSXA.headerTow  = imuTable(:,2);
    outData.RAWIMUSXA.imuError   = imuTable(:,3);
    outData.RAWIMUSXA.imuType    = imuTable(:,4);
    outData.RAWIMUSXA.week       = imuTable(:,5);
    outData.RAWIMUSXA.tow        = imuTable(:,6);
    outData.RAWIMUSXA.epochs     = navsu.time.gps2epochs(outData.RAWIMUSXA.week,outData.RAWIMUSXA.tow);
    outData.RAWIMUSXA.imuStatus  = imuStatus;
    %                           x              y              z
    outData.RAWIMUSXA.acc       = [imuTable(:,10)  -imuTable(:,9) imuTable(:,8)];
    outData.RAWIMUSXA.gyro      = [imuTable(:,13) -imuTable(:,12) imuTable(:,11)];
    
    % scale factors- pulled from OEM7 reference manual, but I have to
    % exclude the 1e-9 factor for the acceleration to look correct. Not
    % clear about the gyro factor.
    %         outData.RAWIMUSXA.acc = outData.acc*2/100;  % this looks correct.
    %         outData.RAWIMUSXA.gyro = outData.gyro/100/100;  % this is probably wrong still.
else
    % Just initialize an empty output
    % initialize output structure
    outData.RAWIMUSXA = struct('headerWeek',[],'headerTow',[],...
        'imuError',[],'imuType',[],'week',[],'tow',[],'epochs',[],'imuStatus',[],'acc',[],'gyro',[]);
end

if RANGEA && ~isempty(find(contains(textdata,'#RANGEA'), 1))
    % find range data
    tempi = textdata(contains(textdata,'#RANGEA'));
    
    epochsi = nan(size(tempi,1),1);
    nEpochs = length(epochsi);
    
    nDataAdd = 1e5;
    dataFull = nan(nDataAdd,26);
    nData = size(dataFull,1);
    indData = 1;
    
    if strcmp(fileExt,'.GPS')
        formatHeader = '%s%f%f%s%f%f%f%s%f%f%f';%f%f%f%f%f%f%f%f%s';
        nHeaderEl = 11;
    else
        formatHeader = '%s%f%f%s%f%f%f%s%f%f';%f%f%f%f%f%f%f%f%s';
        nHeaderEl = 10;
    end
    
    
    % Parse each entry
    for idx = 1:nEpochs
        % read the header
        txti = tempi{idx};
        
        parseHeader = textscan(txti(10:end),formatHeader,1,'Delimiter',',;*');
        
        nObs = parseHeader{end};
        
        
        formatFull = [formatHeader repmat('%f%f%f%f%f%f%f%f%f%s',1,nObs)];
        parseFull = textscan(txti(10:end),formatFull,1,'Delimiter',',;*');
        
        weeki = parseFull{5};
        towi  = parseFull{6};
        % Put it in a nice matrix
        
        epochsi(idx) = navsu.time.gps2epochs(weeki,towi);
        
        if nObs == 0
            continue;
        end
        
        obsMat = reshape(parseFull((nHeaderEl+1):end)',10,nObs)';
        
        % Channel tracking status
        statusBin = fliplr(dec2bin(hex2dec([obsMat{:,end}]),32));
        
        trackingState = bin2dec(statusBin(:,5:-1:1));
        svChannelNumber = bin2dec(statusBin(:,10:-1:6));
        phaseLockFlag = bin2dec(statusBin(:,11));
        parityKnownFlag = bin2dec(statusBin(:,12));
        codeLockedFlag = bin2dec(statusBin(:,13));
        correlatorType = bin2dec(statusBin(:,16:-1:14));
        satSystem      = bin2dec(statusBin(:,19:-1:17));
        grouping       = bin2dec(statusBin(:,21));
        signalType     = bin2dec(statusBin(:,26:-1:22));
        primaryL1Channel = bin2dec(statusBin(:,28));
        carrierPhaseMeas = bin2dec(statusBin(:,29));
        digFiltering     = bin2dec(statusBin(:,30));
        prnLock          = bin2dec(statusBin(:,31));
        chanAssignment   = bin2dec(statusBin(:,32));
        
        obsMati = cell2mat(obsMat(1:end,1:9));
        
        % check if we have enough space in our data matrix
        % read the observation data
        if indData+nObs-1 > nData
            % need to add more space
            data2 =nan(nDataAdd+nData,26);
            %                  = size(datai,1);
            data2(1:size(dataFull,1),:) = dataFull;
            dataFull = data2;
            
            nData = size(dataFull,1);
            clear data2;
        end
        
        % add each observation to our data
        dataFull(indData:(indData+nObs-1),:) = [idx*ones(nObs,1) weeki*ones(nObs,1) ...
            towi*ones(nObs,1) obsMati trackingState svChannelNumber ...
            phaseLockFlag parityKnownFlag codeLockedFlag correlatorType ...
            satSystem grouping signalType primaryL1Channel carrierPhaseMeas ...
            digFiltering prnLock chanAssignment];
        
        indData = indData+nObs;
    end
    
    % Pull off extra data
    dataFull = dataFull(1:(indData-1),:);
    epochsFull = navsu.time.gps2epochs(dataFull(:,2),dataFull(:,3));
    epochs = unique(epochsFull);
    epochs = epochsi;
    
    prnConstSig = unique(dataFull(:,[4 19 21]),'rows');
    nPrnConstSig = size(prnConstSig,1);
    
    % Map from NovAtel constellation to SU
    constMap = [0 1 2 3 4 5;   % novatel
        1 2 0 3 0 0];  % su
    
    constUnique = constMap(2,ismember(constMap(1,:),unique(prnConstSig(:,2))));
    GPS_flag = ismember(1,constUnique);
    GLO_flag = ismember(2,constUnique);
    GAL_flag = ismember(3,constUnique);
    
    if isempty(constellations)
        constellations = navsu.readfiles.initConstellation(GPS_flag,GLO_flag,GAL_flag,0,0,0);
    end
    
    prns      = constellations.PRN;
    constInds = constellations.constInds;
    sInds     = constellations.indexes;
    
    outData.RANGEA.obsData = [];
    outData.RANGEA.epochs = epochs;
    % go through each one and add it
    for idx = 1:nPrnConstSig
        prni = prnConstSig(idx,1);
        consti = constMap(2,constMap(1,:) == prnConstSig(idx,2));
        sigNovi = prnConstSig(idx,3);
        
        switch consti
            case 1 % GPS
                sigInds = [0 5 9 14 17 16]';
                sigMap = {'C1C','L1C','D1C','S1C'; % 0 L1CA
                    'C1L','L1L','D1L','S1L'; % 5 L1C(P)
                    'C2S','L2S','D2S','S2S'; % 9 L2C(M)
                    'C2P','L2P','D2P','S2P'; % 14 L2P
                    'C2W','L2W','D2W','S2W'; % 16 L2P(Y)
                    'C5Q','L5Q','D5Q','S5Q'}; % 17 L5(Q)
                
                sigInds = [0 5 17 14 9 16]';
                sigMap = {'C1C','L1C','D1C','S1C'; % 0 L1CA
                    'C1L','L1L','D1L','S1L'; % 5 L1C(P)
                    'C2S','L2S','D2S','S2S'; % 17  L2P(Y)
                    'C2P','L2P','D2P','S2P'; % 14 L2P
                    'C2W','L2W','D2W','S2W'; % 9  L2C(M)
                    'C5Q','L5Q','D5Q','S5Q'}; % 16 L5(Q)
            case 2 % GLO
                sigInds = [0 1 5 6]';
                sigMap = {'C1C','L1C','D1C','S1C';  % 0 L1CA
                    'C2C','L2C','D2C','S2C';  % 1 L2CA
                    'C2P','L2P','D2P','S2P';  % 5 L2P
                    'C3Q','L3Q','D3Q','S3Q';};% 6 L3(Q)
                prni = prni-37;
            case 3 % GAL
                sigInds = [2 12 17 20 7]';
                sigMap = {'C1C','L1C','D1C','S1C'; % 2 E1C
                    'C5Q','L5Q','D5Q','S5Q'; % 12 E5A(Q)
                    'C7Q','L7Q','D7Q','S7Q'; % 17 E5B(Q)
                    'C8Q','L8Q','D8Q','S8Q'; % 20 E5ALTBOC(Q)
                    'C6C','L6C','D6C','S6C'};% 7 E6C
        end
        
        % check if we have this
        if ~ismember([prni consti],[prns' constInds'],'rows')
            continue;
        end
        
        sIndi = sInds(prns == prni & constInds == consti);
        
        sigNamesi = sigMap(sigInds == sigNovi,:);
        
        if ~isfield(outData.RANGEA.obsData,sigNamesi{1})
            % make all of them
            for jdx = 1:length(sigNamesi)
                outData.RANGEA.obsData.(sigNamesi{jdx}) = nan(length(prns),nEpochs);
                outData.RANGEA.tLock.(sigNamesi{jdx})   = nan(length(prns),nEpochs);
                outData.RANGEA.carrPhaseHalfCycle.(sigNamesi{jdx})   = zeros(length(prns),nEpochs);
            end
        end
        
        % map to the right things - code, carrier, snr measurements
        indsi = find(dataFull(:,4) == prnConstSig(idx,1) &  ...
            dataFull(:,19) == prnConstSig(idx,2) & ...
            dataFull(:,21) == prnConstSig(idx,3));
        datai = dataFull(indsi ,:);
        epochsi = epochsFull(indsi);
        
        % Everything in its right place
        [~,indsLeft] =ismember(epochsi,epochs);
        
        % code phase
        outData.RANGEA.obsData.(sigNamesi{1})(sIndi,indsLeft) = datai(:,6);
        % carrier phase
        outData.RANGEA.obsData.(sigNamesi{2})(sIndi,indsLeft) = -datai(:,8);
        % doppler
        outData.RANGEA.obsData.(sigNamesi{3})(sIndi,indsLeft) = datai(:,10);
        % SNR
        outData.RANGEA.obsData.(sigNamesi{4})(sIndi,indsLeft) = datai(:,11);
        % Save carrier lock time as well
        outData.RANGEA.tLock.(sigNamesi{2})(sIndi,indsLeft)   = datai(:,12);
        % Tracking status
        outData.RANGEA.carrPhaseHalfCycle.(sigNamesi{2})(sIndi,indsLeft)   = datai(:,23);
    end
    
    
    if correctHalfCycle
        % Adjust all ranges for the half cycle offset applied by the
        % receiver for when it incorrectly initialized the phase
        fields = fieldnames(outData.RANGEA.obsData);
        
        for idx = 1:length(fields)
            outData.RANGEA.obsData.(fields{idx}) =  outData.RANGEA.obsData.(fields{idx})+...
                outData.RANGEA.carrPhaseHalfCycle.(fields{idx})*0.5;
        end
        
    end
else
    outData.RANGEA = struct('obsData',[],'epochs',[],'tLock',[],'carrPhaseHalfCycle',[]);
end

if RAWIMUA && ~isempty(find(contains(textdata,'#RAWIMUA'), 1))
    % find imudata
    formati = '%s%s%f%f%s%f%f%f%s%f%f%f%f%f%f%f%f%f%f%s';
    
    tempi = textdata(contains(textdata,'#RAWIMUA'));
    tempi = [tempi{:}];
    
    imuTable = textscan(tempi(2:end),formati,'Delimiter',',;*#');
    
    outData.RAWIMUA.headerWeek = imuTable{6};%[imuTable{:,2}]';
    outData.RAWIMUA.headerTow  = imuTable{7};%[imuTable{:,3}]';
    outData.RAWIMUA.week       = imuTable{11};%[imuTable{:,10}]';
    outData.RAWIMUA.tow        = imuTable{12};%[imuTable{:,11}]';
    outData.RAWIMUA.epochs     = gps2epochs(outData.RAWIMUA.week,outData.RAWIMUA.tow);
    outData.RAWIMUA.imuStatus  = imuTable{13};%imuTable(:,12);
    %                                     x              y              z
    outData.RAWIMUA.acc       = [imuTable{16}  -imuTable{15} imuTable{14}];
    outData.RAWIMUA.gyro      = [imuTable{19}  -imuTable{18} imuTable{17}];
    
else
    % initialize empty output
    outData.RAWIMUA = struct('headerWeek',[],'headerTow',[],...
        'week',[],'tow',[],'epochs',[],'imuStatus',[],'acc',[],'gyro',[]);
end

% end





end