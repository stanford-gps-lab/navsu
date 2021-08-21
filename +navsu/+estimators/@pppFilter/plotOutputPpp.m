function plotOutputPpp(obj,outputs,varargin)
%
p = inputParser;

p.addParameter('truePosEcef',[]);
p.addParameter('truthFile',[]);

% parse the results
parse(p, varargin{:});
res = p.Results;
truthFile            = res.truthFile;      %

% Plot the position and clock bias in ENU
%% Pull out saved values
xyz = [outputs.pos]';
xyzCov = cat(3,outputs.covPos);
epochs = [outputs.epoch]';
b = zeros(size(epochs));
[week,tow] = navsu.time.epochs2gps(epochs);

epochs0 = epochs;
tinds0  = 1:length(epochs);

llh0 = navsu.geo.xyz2llh(xyz(1,:));

utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));

% Convert estimated position to enu
enu    = nan(size(xyz));
enuStd = nan(size(xyz));
for idx = 1:size(xyz,1)
    [enu(idx,1),enu(idx,2),enu(idx,3)] = ...
        navsu.thirdparty.cart2utm(xyz(idx,1),xyz(idx,2),xyz(idx,3),utmZone);
    
    [~,RxyzEnu] = navsu.geo.xyz2enu([0 0 0],llh0(1)*pi/180,llh0(2)*pi/180);
    
    covEnui = RxyzEnu'*squeeze(xyzCov(:,:,idx))*RxyzEnu;
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


%% plot position error
if ~isempty(truth)
    figure;
    ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]);
    
    yplot = [estim.enuPos-truth.enuPosInterp b*navsu.constants.c]';
    yplotStd = [estim.enuStd b*navsu.constants.c]';
    tplot = (epochs-epochs(1))/60; % minutes
    ylabels = {'East [m]' 'North [m]' 'Up [m]' 'Clock [m]'};
    for idx = 1:3
        axes(ha(idx))
        
        plot(tplot,abs(yplot(idx,:)))
        hold on;
        plot(tplot,2*yplotStd(idx,:),'k')
%         plot(tplot,3*-yplotStd(idx,:),'k')
        
        ylabel(ylabels{idx})
        grid on
        if idx < 3
            xticklabels('')
        else
            xlabel('Time [min]')
        end
    end
end

%%
if ~isempty(truth)
    % Plot the ground track and altitude
    figure;
    subplot(4,1,1:3);
   
    % Median position of the ground track
    posMedEnu = nanmedian(estim.enuPos(:,1:2));
    
    enuPlotEst = (estim.enuPos(:,1:2)-posMedEnu);
    enuPlotTrue = (truth.enuPosInterp(:,1:2)-posMedEnu);
    
    s = plot([enuPlotEst(:,1) enuPlotTrue(:,1)],[enuPlotEst(:,2) enuPlotTrue(:,2)],'markersize',10,'linewidth',2);
    
    linestyles = {'none','none'};
    markers = {'.','o'};
    markersizes = [10 6];
    linewidths = [1 2];
    % Add time to data tips
    for idx = 1:length(s)
        row = dataTipTextRow('tIndex',tinds0);
        s(idx).DataTipTemplate.DataTipRows(3) = row;
        s(idx).Marker = markers{idx};
        s(idx).LineStyle = linestyles{idx};
        s(idx).MarkerSize = markersizes(idx);
        s(idx).LineWidth  = linewidths(idx);
    end
    grid on
    legend('Estimated position','Truth position','location','best')
    xlabel('East position [m]')
    ylabel('North position [m]')
    
    subplot(4,1,4);
    s = plot(tinds0,[estim.enuPos(:,3) truth.enuPosInterp(:,3)]);
    
    % Add time to data tips
    for idx = 1:length(s)
        row = dataTipTextRow('tIndex',tinds0);
        s(idx).DataTipTemplate.DataTipRows(3) = row;
    end
    grid on;
    xlabel('Time index')
    ylabel('Altitude [m]')
    
end





%% plot GNSS residuals
% Pull out only the GNSS residuals
residsData = [outputs.resids]';
residsFull = cat(1,residsData.resids);
measIdFull = cat(1,residsData.measId);
epochsFull = cat(1,residsData.epochs);

residsType = cat(1,measIdFull.TypeID);


%%
indsGnss = find(residsType == navsu.internal.MeasEnum.GNSS);

residsGnss = residsFull(indsGnss);
measIdGnss = measIdFull(indsGnss);
epochsGnss = epochsFull(indsGnss);

prnResids = cat(1,measIdGnss.prn);
constResids = cat(1,measIdGnss.const);
freqResids  = cat(1,measIdGnss.freq);
subtypeResids = cat(1,measIdGnss.subtype);

% Individual measurement types
measIdMat = [prnResids constResids freqResids uint8(subtypeResids)];
measIdUn = unique(measIdMat,'rows');

epochMin = min(epochsGnss);

figure;
ha = navsu.thirdparty.tightSubplot(3,1,0.05,[0.1 0.1],[0.07 0.05]);
axes(ha(1))

% Code phase residuals
indsUnPr = measIdUn(measIdUn(:,4) == uint8(navsu.internal.MeasEnum.Code),:);
nPr = size(indsUnPr,1);
epochsUn = unique(epochsGnss);
nEpochs = length(epochsUn);

prResidsPlot = nan(nPr,nEpochs);

for idx = 1:size(indsUnPr)
    indsi = find(ismember(measIdMat,indsUnPr(idx,:),'rows'));
    
    epochsi = epochsGnss(indsi);
    prResidsi     = residsGnss(indsi);
    
    [~,ixb] = ismember(epochsi,epochsUn);
    
    prResidsPlot(idx,ixb) = prResidsi;
end
tPlot = (epochsUn-epochsUn(1))/60;
s = plot(tPlot,prResidsPlot,'.');
for idx = 1:length(s)
    s(idx).DataTipTemplate.DataTipRows(1).Label = 't';
    s(idx).DataTipTemplate.DataTipRows(2).Label = 'resid';
    row = dataTipTextRow('PRN',double(indsUnPr(idx,1)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(3) = row;
    row = dataTipTextRow('const',double(indsUnPr(idx,2)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(4) = row;
    row = dataTipTextRow('sig',double(indsUnPr(idx,3)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(5) = row;
end

xlim([0 (max(epochsGnss)-min(epochsGnss))/60]); grid on;
ylabel('Code phase residuals [m]')
title('Measurement residuals over time')

% Plot carrier phase residuals
axes(ha(2))
% Carrier phase residuals
indsUnPh = measIdUn(measIdUn(:,4) == uint8(navsu.internal.MeasEnum.Carrier),:);
nPr = size(indsUnPh,1);
epochsUn = unique(epochsGnss);
nEpochs = length(epochsUn);

phResidsPlot = nan(nPr,nEpochs);

for idx = 1:size(indsUnPh,1)
    indsi = find(ismember(measIdMat,indsUnPh(idx,:),'rows'));
    
    epochsi = epochsGnss(indsi);
    prResidsi     = residsGnss(indsi);
    
    [~,ixb] = ismember(epochsi,epochsUn);
    
    phResidsPlot(idx,ixb) = prResidsi;
end
tPlot = (epochsUn-epochsUn(1))/60;
s = plot(tPlot,phResidsPlot);
for idx = 1:length(s)
    s(idx).DataTipTemplate.DataTipRows(1).Label = 't';
    s(idx).DataTipTemplate.DataTipRows(2).Label = 'resid';
    row = dataTipTextRow('PRN',double(indsUnPh(idx,1)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(3) = row;
    row = dataTipTextRow('const',double(indsUnPh(idx,2)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(4) = row;
    row = dataTipTextRow('sig',double(indsUnPh(idx,3)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(5) = row;
end

xlim([0 (max(epochsGnss)-min(epochsGnss))/60]); grid on;
ylabel('Carrier phase residuals [m]')
ylim([-0.1 0.1])

% Plot doppler residuals
axes(ha(3))
% Doppler residuals
indsUnDopp = measIdUn(measIdUn(:,4) == uint8(navsu.internal.MeasEnum.Doppler),:);
nPr = size(indsUnDopp,1);
epochsUn = unique(epochsGnss);
nEpochs = length(epochsUn);

doppResidsPlot = nan(nPr,nEpochs);

for idx = 1:size(indsUnDopp,1)
    indsi = find(ismember(measIdMat,indsUnDopp(idx,:),'rows'));
    
    epochsi = epochsGnss(indsi);
    prResidsi     = residsGnss(indsi);
    
    [~,ixb] = ismember(epochsi,epochsUn);
    
    doppResidsPlot(idx,ixb) = prResidsi;
end
tPlot = (epochsUn-epochsUn(1))/60;
s = plot(tPlot,doppResidsPlot,'.');
for idx = 1:length(s)
    s(idx).DataTipTemplate.DataTipRows(1).Label = 't';
    s(idx).DataTipTemplate.DataTipRows(2).Label = 'resid';
    row = dataTipTextRow('PRN',double(indsUnDopp(idx,1)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(3) = row;
    row = dataTipTextRow('const',double(indsUnDopp(idx,2)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(4) = row;
    row = dataTipTextRow('sig',double(indsUnDopp(idx,3)).*ones(nEpochs,1));
    s(idx).DataTipTemplate.DataTipRows(5) = row;
end

xlim([0 (max(epochsGnss)-min(epochsGnss))/60]); grid on;
xlabel('Minutes into run')
ylabel('Doppler residuals [m]')


% Also plot summaries of the residuals
prnConstList = sortrows(unique(measIdMat(:,1:2),'rows'),2);

figure;
ha = navsu.thirdparty.tightSubplot(3,1,0.05,[0.1 0.1],[0.07 0.05]);
axes(ha(1))
% pseudorange means
means = nanmean(prResidsPlot,2);
stds  = nanstd(prResidsPlot');
[~,xPlotPr] = ismember(indsUnPr(:,1:2),prnConstList,'rows');
s = errorbar(xPlotPr,means,stds,'o');
grid on;
row = dataTipTextRow('PRN',double(indsUnPr(:,1)));
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('const',double(indsUnPr(:,2)));
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('sig',double(indsUnPr(:,3)));
s.DataTipTemplate.DataTipRows(end+1) = row;
ylabel('Code resids [m]')


axes(ha(2))
% carrier phase means
means = nanmean(phResidsPlot,2);
stds  = nanstd(phResidsPlot');
[~,xPlotPh] = ismember(indsUnPh(:,1:2),prnConstList,'rows');
s = errorbar(xPlotPh,means,stds,'o');
grid on;
row = dataTipTextRow('PRN',double(indsUnPh(:,1)));
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('const',double(indsUnPh(:,2)));
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('sig',double(indsUnPh(:,3)));
s.DataTipTemplate.DataTipRows(end+1) = row;
ylabel('Carrier residuals [m]')


axes(ha(3))
% doppler means
means = nanmean(doppResidsPlot,2);
stds  = nanstd(doppResidsPlot');
[~,xPlotDopp] = ismember(indsUnDopp(:,1:2),prnConstList,'rows');
s = errorbar(xPlotDopp,means,stds,'o');
grid on;
row = dataTipTextRow('PRN',double(indsUnDopp(:,1)));
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('const',double(indsUnDopp(:,2)));
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('sig',double(indsUnDopp(:,3)));
s.DataTipTemplate.DataTipRows(end+1) = row;
ylabel('Doppler residuals [m]')
xlabel('Minutes into run')


%% plot measurements that were removed?
if 1
    measRemovedStruc = cat(1,outputs.measRemoved);
    measRemId = cat(1,measRemovedStruc.id);
    % only keeping GNSS measurements
    indsGnss = find(cat(1,measRemId.TypeID) == navsu.internal.MeasEnum.GNSS);
    
    measRemReas = cat(1,measRemovedStruc.reason);
    measRemEpoch = cat(1,measRemovedStruc.epoch);
    measRemReas= measRemReas(indsGnss);
    measRemEpoch = measRemEpoch(indsGnss);
    measRemId = measRemId(indsGnss);
    
    measRemSubtype = cat(1,measRemId.subtype);
    
    measRemPrnConst = [cat(1,measRemId.prn) cat(1,measRemId.const)];
    
    prns = outputs(1).residsInfo.rangeInfo.PRN(1,:)';
    constInds = outputs(1).residsInfo.rangeInfo.constInds(1,:)';
    sInds = (1:length(prns))';
    epochs = [outputs.epoch]';
    
    [~,yInds] = ismember(measRemPrnConst(:,1:2),[prns constInds],'rows');
    [~,xInds] =  ismember(measRemEpoch,epochs);
    
    yPlotMat = repmat(sInds,1,length(epochs));
    xPlotMat = repmat(1:length(epochs),length(sInds),1);
    
    figure; hold on;
    
    markers = {'g.','c.','ro','r.','b^','k'};
    markerSize = [5 2 5 8 2 2];
    legText = {'Meas Used','Low Elevation','Code Removed','Carr Removed','Cycle Slip','Number of satellites used'};
    for idx = 1:length(legText)
        plot(-10,-10,markers{idx},'markerSize',markerSize(idx))
    end
    
    % general availability of measurements
    prnConstEpochs = unique([cat(1,measIdGnss.prn) cat(1,measIdGnss.const) epochsGnss],'rows');
    [~,yIndsAvail] = ismember(prnConstEpochs(:,[1 2]),[prns constInds],'rows');
    [~,xIndsAvail] = ismember(prnConstEpochs(:,3) ,epochs);
    
    for idx = 1:length(legText)
        switch idx
            case 1 % Meas used
                linInd = sub2ind(size(yPlotMat),yIndsAvail,xIndsAvail);
                
            case 2 % low elevation
                indsi = find(measRemReas == 1);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
                
            case 3 % code removed
                indsi = find(measRemReas == 2 & measRemSubtype == navsu.internal.MeasEnum.Code);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
                
            case 4 % carrier removed
                indsi = find(measRemReas == 2 & measRemSubtype == navsu.internal.MeasEnum.Carrier);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
            case 5
                % cycle slip
                indsi = find( measRemReas == 3);
                linInd = sub2ind(size(yPlotMat),yInds(indsi),xInds(indsi));
        end
        
        xPloti = nan(size(xPlotMat));
        xPloti(linInd) = xPlotMat(linInd);
        yPloti = nan(size(yPlotMat));
        yPloti(linInd) = yPlotMat(linInd);
        
        plot(xPloti',yPloti',markers{idx},'markerSize',markerSize(idx));
    end
    %     plot(xPlotMat(1,:),nSatsUsed,'k');
    xlim([min(min(xPlotMat))-0.5 max(max(xPlotMat))+0.5])
    ylim([min(min(yPlotMat))-0.5 max(max(yPlotMat))+0.5])
    legend(legText);
    
end


%% Plot sky plot
% collect azimuths and elevations
residsData = [outputs.resids]';
elFull = cat(1,residsData.el);
azFull = cat(1,residsData.az);

prnConstFull = cat(1,residsData.prnConstInds);
epochsFull = cat(1,residsData.epochsElAz);

satsUn = unique(prnConstFull,'rows');
epochsUn = unique(epochsFull);

el = nan(length(satsUn),length(epochsUn));
az = nan(length(satsUn),length(epochsUn));

for idx = 1:size(satsUn,1)
   indsi =  find(prnConstFull(:,1) == satsUn(idx,1) & prnConstFull(:,2) == satsUn(idx,2));
    
   eli = elFull(indsi);
   azi = azFull(indsi);
   epochsi = epochsFull(indsi);
   
   [~,ixb] = ismember(epochsi,epochsUn);
   el(idx,ixb) = eli;
   az(idx,ixb) = azi;
end

figure
navsu.geo.skyPlot(az'*180/pi,el'*180/pi,'label',true,'prn',satsUn(:,1),...
    'const',satsUn(:,2),'linestyle','o')


%%
% plot flex states
flex = cat(1,outputs.flexStates);

flexStates = cat(1,flex.states);
flexInfo   = cat(1,flex.info);
flexStd    = cat(1,flex.std);
flexEpochs = cat(1,flex.epochs);

%% plot ambiguities
indsAmb = find(flexInfo(:,3) == 1);

% Collect all of the ambiguities to plot
flexInfoAmb = flexInfo(indsAmb,:);
flexStatesAmb = flexStates(indsAmb,:);
flexStdAmb    = flexStd(indsAmb);
flexEpochs    = flexEpochs(indsAmb);

ambsUn = unique(flexInfoAmb,'rows');
epochsUn = unique(flexEpochs);

ambStatePlot = nan(length(ambsUn),length(epochsUn));
ambStdPlot   = nan(length(ambsUn),length(epochsUn));

for idx = 1:length(ambsUn)
    indsi = find(ismember(flexInfoAmb,ambsUn(idx,:),'rows'));
    
    statesi = flexStatesAmb(indsi);
    stdi    = flexStdAmb(indsi);
    epochsi = flexEpochs(indsi);
    
    [~,ixb] = ismember(epochsi,epochsUn);
    
    ambStatePlot(idx,ixb) = statesi;
    ambStdPlot(idx,ixb) = stdi;
end

colors = lines(length(ambsUn));

figure
s = plot(ambStatePlot','linewidth',2);
hold on
for idx = 1:length(s)
   % add data tip information 
    s(idx).DataTipTemplate.DataTipRows(1).Label = 't';
    s(idx).DataTipTemplate.DataTipRows(2).Label = 'resid';
    row = dataTipTextRow('PRN',double(ambsUn(idx,1)).*ones(length(epochsUn),1));
    s(idx).DataTipTemplate.DataTipRows(3) = row;
    row = dataTipTextRow('const',double(ambsUn(idx,2)).*ones(length(epochsUn),1));
    s(idx).DataTipTemplate.DataTipRows(4) = row;
     row = dataTipTextRow('sig',double(ambsUn(idx,4)).*ones(length(epochsUn),1));
    s(idx).DataTipTemplate.DataTipRows(4) = row;
    s(idx).Color = colors(idx,:);
end
ylabel('Integer ambiguities')
% 
% % Add standard dev lol
% s = plot(ambStatePlot'+ambStdPlot');
% for idx = 1:length(s)
%    % add data tip information 
%     s(idx).DataTipTemplate.DataTipRows(1).Label = 't';
%     s(idx).DataTipTemplate.DataTipRows(2).Label = 'std+resid';
%     row = dataTipTextRow('PRN',double(ambsUn(idx,1)).*ones(length(epochsUn),1));
%     s(idx).DataTipTemplate.DataTipRows(3) = row;
%     row = dataTipTextRow('const',double(ambsUn(idx,2)).*ones(length(epochsUn),1));
%     s(idx).DataTipTemplate.DataTipRows(4) = row;
%      row = dataTipTextRow('sig',double(ambsUn(idx,4)).*ones(length(epochsUn),1));
%     s(idx).DataTipTemplate.DataTipRows(4) = row;
%     s(idx).Color = colors(idx,:);
% end
% s = plot(ambStatePlot'-ambStdPlot');
% for idx = 1:length(s)
%    % add data tip information 
%     s(idx).DataTipTemplate.DataTipRows(1).Label = 't';
%     s(idx).DataTipTemplate.DataTipRows(2).Label = 'std+resid';
%     row = dataTipTextRow('PRN',double(ambsUn(idx,1)).*ones(length(epochsUn),1));
%     s(idx).DataTipTemplate.DataTipRows(3) = row;
%     row = dataTipTextRow('const',double(ambsUn(idx,2)).*ones(length(epochsUn),1));
%     s(idx).DataTipTemplate.DataTipRows(4) = row;
%      row = dataTipTextRow('sig',double(ambsUn(idx,4)).*ones(length(epochsUn),1));
%     s(idx).DataTipTemplate.DataTipRows(4) = row;
%     s(idx).Color = colors(idx,:);
% end

%% Plot the dop 
resids2 = cat(1,outputs.resids);

dop = cat(2,resids2.dop);
epochsDop = cat(2,resids2.dopEpoch);
dop(dop == 0) = nan;

figure;
plot(epochsDop,dop(1:3,:))
legend('E','N','U')
xlabel('epochs')
ylabel('DOP')

end



























