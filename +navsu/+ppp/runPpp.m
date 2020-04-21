function outData = runPpp(filter,obsGnss,corrData,varargin)

%% Sort all of the GNSS and IMU measurements
% 1 = GNSS, 2 = IMU
%                      EPOCH      |        GNSS/IMU            |  INDEX WITHIN TYPE
obsInfo = sortrows([obsGnss.epochs 1*ones(size(obsGnss.epochs)) [1:size(obsGnss.epochs)]'],1);

indStart = min(find(obsInfo(:,2) == 1));
obsInfo = obsInfo(indStart:end,:);

nEpochs = size(obsInfo,1);

% Initialize the waitbar
runTimeStart = tic;
pctDone = 0;
h = waitbar(0,'0 Percent Complete');

%% Data to save
outData = [];

%% Run the loop
for tdx = 1:nEpochs
    % Pull the measurement from the full list
    obsi = navsu.ppp.stripMeas({obsGnss},obsInfo(tdx,2),obsInfo(tdx,3));
    
    epochi = obsInfo(tdx,1);
    
    % GNSS measurements
    if ~any(obsi.range.obs(:))
        updated = false;
        continue;
    end
    
    % if not initialized, try to initialize :)
    if filter.initialized
        % Do the ppp update
        % Do the time and measurement updates
        filter.update(epochi,obsi,corrData);
        
    else
        % need to initialize
        filter.initialize(corrData,'gnssMeas',obsi);
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

