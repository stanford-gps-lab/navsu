function outData = runPpp(filter,obs,corrData,varargin)

%% 

nEpochs = size(obs,1);

% Initialize the waitbar
runTimeStart = tic;
pctDone = 0;
h = waitbar(0,'0 Percent Complete');

%% Data to save
outData = [];

%% Run the loop
for tdx = 1:nEpochs
    % Pull the measurement from the full list
    obsi = obs{tdx};
    
    epochi = obsi{1}.epochs;
    
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

