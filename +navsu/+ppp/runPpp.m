function outData = runPpp(filter,obs,corrData,varargin)

%% Sync all measurements 
% measTypeEnums = enumeration(navsu.internal.MeasEnum.GNSS);
% measStruc = cell2struct(repmat({[]},length(measTypeEnums),1),cellstr(measTypeEnums));

% go through each measurement struct given and put it in here
epochsFull =[];
for idx = 1:length(obs)
    epochsFull = [epochsFull; obs{idx}.epochs];
end
epochs = unique(epochsFull);
obsMap = nan(length(epochs),length(obs));
for idx = 1:length(obs)
    [~,ixb] = ismember(epochs,obs{idx}.epochs);
    obsMap(:,idx) = ixb;
end

%%
nEpochs = length(epochs);

% Initialize the waitbar
runTimeStart = tic;
pctDone = 0;
h = waitbar(0,'0 Percent Complete');

%% Data to save
outData = [];

%% Run the loop
for tdx = 1:nEpochs
    % Pull the measurement from the full list
    obsi = navsu.ppp.stripMeas(tdx,obs,obsMap);
    
    epochi = epochs(tdx);
    
    % if not initialized, try to initialize :)
    if filter.initialized
        % Do the time and measurement updates
        filter.update(epochi,obsi,corrData);
    else
        % need to initialize
        filter.initialize(corrData,obsi);
    end
    
    % Save some things for output
    outData = filter.saveState(outData,epochi,obsi);
    
    % Update the waitbar
    if mod(floor(tdx/nEpochs*100),1) == 0 && floor(tdx/nEpochs*100) > pctDone
        tElapsed = toc(runTimeStart);
        tRemaining = tElapsed*(nEpochs-tdx)./tdx;
        pctDone =  floor(tdx/nEpochs*100);
        waitbar(pctDone/100,h,[num2str(pctDone) '% Complete, ' num2str(tRemaining/60,'%5.2f') ' Minutes Remaining']);
    end
end

close(h);


end

