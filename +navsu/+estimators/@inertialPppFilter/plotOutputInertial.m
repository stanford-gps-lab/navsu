function plotOutputInertial(outputs,varargin)
%
p = inputParser;

p.addParameter('truePosEcef',[]);
p.addParameter('truthFile',[]);

% parse the results
parse(p, varargin{:});
res = p.Results;
truePosEcef          = res.truePosEcef;          % truth position to compare to
truthFile            = res.truthFile;      %

% Plot the position and clock bias in ENU
%% Pull out saved values
xyz = [outputs.pos]';
xyzCov = cat(3,outputs.covPos);
epochs = [outputs.epoch]';
b = zeros(size(epochs));
Rbe = cat(3,outputs.R_b_e);


llh0 = navsu.geo.xyz2llh(xyz(1,:));

utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));

% Convert estimated position and attitude to ENU
enu    = nan(size(xyz));
enuStd = nan(size(xyz));
attEnu = nan(size(xyz));
for idx = 1:size(xyz,1)
    [enu(idx,1),enu(idx,2),enu(idx,3)] = ...
        navsu.thirdparty.cart2utm(xyz(idx,1),xyz(idx,2),xyz(idx,3),utmZone);
    
    [~,RxyzEnu] = navsu.geo.xyz2enu([0 0 0],llh0(1)*pi/180,llh0(2)*pi/180);
    
    covEnui = RxyzEnu'*squeeze(xyzCov(:,:,idx))*RxyzEnu;
    enuStd(idx,:) = sqrt(diag(covEnui));
    
    Rbenu = RxyzEnu*Rbe(:,:,idx);
    
    attEnu(idx,:) = navsu.geo.dcm2euler123(Rbenu);
    
end
estim.pos = xyz;
estim.epochs = epochs;
estim.enuPos = enu;
estim.enuStd = enuStd;

if ~isempty(truthFile)
    % Parse the truth data
    [~,posTruth, epochsTruth,cov,~,~,velEnu,stdEnu] = navsu.internal.parsePosSolFile(truthFile);
    
    truth.pos = posTruth;
    truth.epochs = epochsTruth;
    truth.enuPos = nan(size(posTruth));
    truth.enuStd = stdEnu;
    
    llh0 = navsu.geo.xyz2llh(posTruth(1,:));
    % Convert truth data to ENU
    utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));
    for idx = 1:size(posTruth,1)
        [truth.enuPos(idx,1),truth.enuPos(idx,2),truth.enuPos(idx,3)] = ...
            navsu.thirdparty.cart2utm(posTruth(idx,1),posTruth(idx,2),posTruth(idx,3),utmZone);
    end
    
    % Interpolate truth data to match the estimated data :)
    truthEnuInterp = nan(size(estim.enuPos));
    truthXyzInterp = nan(size(estim.enuPos));
    if isempty(truth.epochs)
        truthEnuInterp = repmat(truth.enuPos,size(estim.enuPos,1),1);
        truthXyzInterp = repmat(truth.pos,size(estim.enuPos,1),1);
    else
        % Pad truth gaps with NaNs
        indsGap = find(diff(truth.epochs) > 1);
        truth.enuPos(indsGap,:) = nan(length(indsGap),3);
        truth.pos(indsGap,:) = nan(length(indsGap),3);
        
        interpMethod = 'linear';
        truthEnuInterp(:,1) = interp1(truth.epochs,truth.enuPos(:,1),estim.epochs,interpMethod);
        truthEnuInterp(:,2) = interp1(truth.epochs,truth.enuPos(:,2),estim.epochs,interpMethod);
        truthEnuInterp(:,3) = interp1(truth.epochs,truth.enuPos(:,3),estim.epochs,interpMethod);
        
        truthXyzInterp(:,1) = interp1(truth.epochs,truth.pos(:,1),estim.epochs,interpMethod);
        truthXyzInterp(:,2) = interp1(truth.epochs,truth.pos(:,2),estim.epochs,interpMethod);
        truthXyzInterp(:,3) = interp1(truth.epochs,truth.pos(:,3),estim.epochs,interpMethod);
    end
    
    truth.enuPosInterp = truthEnuInterp;
    truth.xyzPosInterp = truthXyzInterp;
else
    % There is no truth.  But this you must learn for yourself.
    truth = [];
end


%% plot

if ~isempty(truth)
    figure;
    ha = navsu.thirdparty.tightSubplot(4,1,0.02,[0.1 0.1],[0.07 0.05]);
    
    yplot = [estim.enuPos-truth.enuPosInterp b*navsu.constants.c]';
    yplotStd = [estim.enuStd b*navsu.constants.c]';
    tplot = (epochs-epochs(1))/60; % minutes
    ylabels = {'East [m]' 'North [m]' 'Up [m]' 'Clock [m]'};
    for idx = 1:4
        axes(ha(idx))
        
        plot(tplot,yplot(idx,:))
        hold on;
        plot(tplot,2*yplotStd(idx,:),'k')
        plot(tplot,2*-yplotStd(idx,:),'k')
        
        ylabel(ylabels{idx})
        grid on
        if idx < 4
            xticklabels('')
        else
            xlabel('Time [min]')
        end
    end
    
end



%% plot attitude

figure;
plot(attEnu*180/pi)

end



























