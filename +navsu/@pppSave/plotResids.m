function plotResids(obj)

% Plot residuals over time
figi = obj.figResids;
if ~isvalid(figi)
    obj.figResids = figure;
    figi = obj.figResids;
end
figi.Visible = 'on';
clf;
figi.Position = [200 200 800 630];

indsPr = find(obj.gnssData.range.ind(:,1) == 1);
indsPh = find(obj.gnssData.range.ind(:,1) == 2);

prns = obj.gnssData.PRN(1,:)';

tPlot = [obj.gnssData.epochs-min(obj.gnssData.epochs)];

ha = navsu.thirdparty.tightSubplot(3,1,0.05,[0.1 0.1],[0.07 0.05]);

residsPr = obj.gnssData.range.resids(indsPr,:,:);
residsPr = reshape(residsPr,size(residsPr,1)*size(residsPr,2),size(residsPr,3));

residsPh = obj.gnssData.range.resids(indsPh,:,:);
residsPh = reshape(residsPh,size(residsPh,1)*size(residsPh,2),size(residsPh,3));
constsPh = obj.gnssData.range.constInds(indsPh,:);
constsPh = reshape(constsPh,size(constsPh,1)*size(constsPh,2),1);
% Plot code phase residuals
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
residsDop = obj.gnssData.doppler.resids;
residsDop = reshape(residsDop,size(residsDop,1)*size(residsDop,2),size(residsDop,3));

% plot(tPlot,residsPh(constsPh == 1,:),'-')
plot(tPlot,residsDop)

xlim([0 max(tPlot)]); grid on;
xlabel('Minutes into run')
ylabel('Doppler residuals [m]')
% ylim([-0.1 0.1])


%%
if 0
    figi = obj.figResidsEl;
    if ~isvalid(figi)
        figi = figure;
    end
    figi.Visible = true;
    clf; hold on;
    
    % map elevation matrix to the measurements we have available
    sIndsMap = repmat(1:size(obj.gnssData.rangeResids,2),size(obj.gnssData.rangeResids,1),1);
    sIndsMapPh = sIndsMap(indsPh,:);
    sIndsMapPh = sIndsMapPh(:);
    sIndsMapPr = sIndsMap(indsPr,:);
    sIndsMapPr = sIndsMapPr(:);
    
    elPh = obj.gnssData.el(sIndsMapPh,:)*180/pi;
    elPr = obj.gnssData.el(sIndsMapPr,:)*180/pi;
    
    indsAvailPh = find(any(~isnan(residsPh)'));
    
    indPloti = 1;
    
    plot(elPh(indsAvailPh(indPloti),:)',residsPh(indsAvailPh(indPloti),:)')
    ylim([-0.1 0.1])
    xlim([0 90])

end





end