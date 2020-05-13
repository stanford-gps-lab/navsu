function outData = runPpp(filter,obs,corrData,varargin)

%% Check what measurements are in here and sort out what's usable
obs = filter.checkMeas(obs);

%% Sync all measurements 
[obsMap,epochs] = navsu.ppp.syncMeas(obs);

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
        filter.initialize(obsi,corrData);
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

