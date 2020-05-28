function outData = runPpp(filter,obs,corrData,varargin)

%% Check what measurements are in here and sort out what's usable
obs = filter.checkMeas(obs);

%% Sync all measurements 
[obsMap,epochs] = navsu.ppp.syncMeas(obs);
nEpochs = length(epochs);

% Initialize the progress bar
wb = navsu.internal.loadingBar(nEpochs);

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
    
    if filter.initialized
        % Save some things for output
        outData = filter.saveState(outData,epochi,obsi);
    end
    
    % Update the progress bar    
    wb.update(tdx);
end

% Close the progress bar
wb.close;

end

