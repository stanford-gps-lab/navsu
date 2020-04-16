function outStruc = runPpp(obsGnss,corrData,PARAMS,varargin)


%% Sort all of the GNSS and IMU measurements
% 1 = GNSS, 2 = IMU
%                      EPOCH      |        GNSS/IMU            |  INDEX WITHIN TYPE
obsInfo = sortrows([obsGnss.epochs 1*ones(size(obsGnss.epochs)) [1:size(obsGnss.epochs)]'],1);

indStart = min(find(obsInfo(:,2) == 1));
obsInfo = obsInfo(indStart:end,:);

nEpochs = size(obsInfo,1);

%% Initialize the object (not the solution)
user = navsu.pppFilter;

% Initialize the waitbar
runTimeStart = tic;
pctDone = 0;
h = waitbar(0,'0 Percent Complete');

neverInitialized = true;

%% Run the loop
for tdx = 1:nEpochs
    % Pull the measurement from the full list
    obsi = navsu.ppp.stripMeas({obsGnss},obsInfo(tdx,2),obsInfo(tdx,3));
    
    epochi = obsInfo(tdx,1);
    
    measType = obsInfo(tdx,2);
    
    updated = true;
    
    % GNSS measurements
    if ~any(obsi.range.obs(:))
        updated = false;
        continue;
    end
    
    % if not initialized, try to initialize :)
    if user.initialized        
        % Do the ppp update
        
        % or just do a least squares solution :)
        user.initialize(PARAMS,corrData,'gnssMeas',obsi);
        
    else
        % need to initialize
        user.initialize(PARAMS,corrData,'gnssMeas',obsi);
        
        if user.initialized && neverInitialized
            % IT WORKED!
            outStruc = navsu.pppSave(obsInfo(:,1),length(user.clockBias),obsGnss,...
                obsGnss.epochs,user.INDS_STATE,obsInfo(:,2),'gnssEpochsOnly',false);
            epochLastPlot = obsInfo(1,1);
            
            outStruc.saveState(user,PARAMS,'tdx', tdx)
        end
    end
    
    if updated
        outStruc.saveState(user,PARAMS,'tdx',tdx);
    end
    
    % Update the waitbar
    if mod(floor(tdx/nEpochs*100),1) == 0 && floor(tdx/nEpochs*100) > pctDone
        tElapsed = toc(runTimeStart);
        tRemaining = tElapsed*(nEpochs-tdx)./tdx;
        pctDone =  floor(tdx/nEpochs*100);
        waitbar(pctDone/100,h,[num2str(pctDone) '% Complete, ' num2str(tRemaining/60,'%5.2f') ' Minutes Remaining']);
    end
end

close(h);
outStruc.closeSaveState


end














