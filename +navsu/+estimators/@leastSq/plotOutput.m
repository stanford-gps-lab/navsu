function plotOutput(outputs,varargin)
% function plotOutput(outputs,varargin)
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

llh0 = navsu.geo.xyz2llh(xyz(1,:));

utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));

% Convert estimated position to enu
enu    = nan(size(xyz));
enuStd = nan(size(xyz));
for idx = 1:size(xyz,1)
    [enu(idx,1),enu(idx,2),enu(idx,3)] = ...
        navsu.thirdparty.cart2utm(xyz(idx,1),xyz(idx,2),xyz(idx,3),utmZone);
    
    [~,RxyzEnu] = navsu.geo.xyz2enu([0 0 0],llh0(1)*pi/180,llh0(2)*pi/180);
    
    covEnui = RxyzEnu*squeeze(xyzCov(:,:,idx))*RxyzEnu';
    enuStd(idx,:) = sqrt(diag(covEnui));
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

%% plot residuals
if 0
    residsData = [outputs.resids];
    
    residsRangeInfo = outputs(1).residsInfo.rangeInfo;
    residsDopplerInfo = outputs(1).residsInfo.dopplerInfo;
    
    residsRange = cat(3,residsData.range);
    residsDoppler = cat(3,residsData.doppler);
    
    tPlot = ([residsData.epoch]-min([residsData.epoch]))/60;
    
    indsPr = find(residsRangeInfo.ind(:,1) == 1);
    indsPh = find(residsRangeInfo.ind(:,1) == 2);
    
    residsPr = residsRange(indsPr,:,:);
    residsPr = reshape(residsPr,size(residsPr,1)*size(residsPr,2),size(residsPr,3));
    
    residsPh = residsRange(indsPh,:,:);
    residsPh = reshape(residsPh,size(residsPh,1)*size(residsPh,2),size(residsPh,3));
    
    figure;
    ha = navsu.thirdparty.tightSubplot(3,1,0.05,[0.1 0.1],[0.07 0.05]);
    axes(ha(1))
    plot(tPlot,residsPr,'.')
    xlim([0 max(tPlot)]); grid on;
    ylabel('Code phase residuals [m]')
    title('Measurement residuals over time')
    
    % Plot carrier phase residuals
    axes(ha(2))
    % plot(tPlot,residsPh(constsPh == 1,:),'-')
    plot(tPlot,residsPh(:,:),'-')
    
    xlim([0 max(tPlot)]); grid on;
    % xlabel('Minutes into run')
    ylabel('Carrier phase residuals [m]')
    ylim([-0.1 0.1])
    
    % Plot doppler residuals
    axes(ha(3))
    residsDoppler = reshape(residsDoppler,size(residsDoppler,1)*size(residsDoppler,2),size(residsDoppler,3));
    
    % plot(tPlot,residsPh(constsPh == 1,:),'-')
    plot(tPlot,residsDoppler,'.')
    
    xlim([0 max(tPlot)]); grid on;
    xlabel('Minutes into run')
    ylabel('Doppler residuals [m]')
    
    %% plot measurements that were removed?
    measRemovedStruc = [outputs.measRemoved];
    measRemoved = cat(1,measRemovedStruc.measRemove);
    epochsRemoved = cat(1,measRemovedStruc.epoch);
    
    prns = outputs(1).residsInfo.rangeInfo.PRN(1,:)';
    constInds = outputs(1).residsInfo.rangeInfo.constInds(1,:)';
    sInds = (1:length(prns))';
    epochs = [residsData.epoch]';
    
    [~,yInds] = ismember(measRemoved(:,1:2),[prns constInds],'rows');
    [~,xInds] =  ismember(epochsRemoved,epochs);
    
    yPlotMat = repmat(sInds,1,length(epochs));
    xPlotMat = repmat(1:length(epochs),length(sInds),1);
    
    figure; hold on;
    
    
    markers = {'g.','c.','ro','r.','b^','k'};
    markerSize = [5 2 5 8 2 2];
    legText = {'Meas Used','Low Elevation','Code Removed','Carr Removed','Cycle Slip','Number of satellites used'};
    for idx = 1:length(legText)
        plot(-10,-10,markers{idx},'markerSize',markerSize(idx))
    end
    rangeResids = cat(3,residsData.range);
    
    for idx = 1:length(legText)
        switch idx
            case 1 % Meas used
                linInd = find(squeeze(sum(~isnan(rangeResids),1)) > 0);
                
                temp = zeros(size(xPlotMat));
                temp(linInd) = 1;
                nSatsUsed = sum(temp,1);
                'fdaf';
                %             nSatsUsed =
            case 2 % low elevation
                indsi = find(measRemoved(:,5) == 1);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
                
            case 3 % code removed
                indsi = find(measRemoved(:,5) == 2 & measRemoved(:,4) == 1);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
                
            case 4 % carrier removed
                indsi = find(measRemoved(:,5) == 2 & measRemoved(:,4) == 2);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
            case 5
                % cycle slip
                indsi = find( measRemoved(:,5) == 3);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
        end
        
        xPloti = nan(size(xPlotMat));
        xPloti(linInd) = xPlotMat(linInd);
        yPloti = nan(size(yPlotMat));
        yPloti(linInd) = yPlotMat(linInd);
        
        plot(xPloti',yPloti',markers{idx},'markerSize',markerSize(idx));
    end
    plot(xPlotMat(1,:),nSatsUsed,'k');
    xlim([min(min(xPlotMat))-0.5 max(max(xPlotMat))+0.5])
    ylim([min(min(yPlotMat))-0.5 max(max(yPlotMat))+0.5])
    legend(legText);
    
end
end

























































