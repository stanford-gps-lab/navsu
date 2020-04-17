function plotOutput(outputs,varargin)
% 
p = inputParser;

p.addParameter('truePosEcef',[]);

% parse the results
parse(p, varargin{:});
res = p.Results;
truePosEcef          = res.truePosEcef;          % truth position to compare to

% Plot the position and clock bias in ENU
xyz = [outputs.pos]';
epochs = [outputs.epoch]';
b = nan(size(epochs));

llh0 = navsu.geo.xyz2llh(xyz(1,:));

utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));

enu = nan(size(xyz));
for idx = 1:size(xyz,1)
    [enu(idx,1),enu(idx,2),enu(idx,3)] = ...
        navsu.thirdparty.cart2utm(xyz(idx,1),xyz(idx,2),xyz(idx,3),utmZone);
end

%
if ~isempty(truePosEcef) 
    truePosEnu = nan(1,3);
    [truePosEnu(1), truePosEnu(2), truePosEnu(3)] = navsu.thirdparty.cart2utm(truePosEcef(1),truePosEcef(2),truePosEcef(3),utmZone);
end


%% plot
figure; 
ha = navsu.thirdparty.tightSubplot(4,1,0.02,[0.1 0.1],[0.07 0.05]);

if ~isempty(truePosEnu)
    compPosEnu = truePosEnu;
else
   compPosEnu = enu(end,:); 
end

yplot = [enu-compPosEnu b*navsu.constants.c]';
tplot = (epochs-epochs(1))/60; % minutes
ylabels = {'East [m]' 'North [m]' 'Up [m]' 'Clock [m]'};
for idx = 1:4
    axes(ha(idx))
    plot(tplot,yplot(idx,:))
    ylabel(ylabels{idx})
    grid on
    if idx < 4
       xticklabels('') 
    else
       xlabel('Time [min]') 
    end
end

%% plot residuals
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
plot(tPlot,residsDoppler)

xlim([0 max(tPlot)]); grid on;
xlabel('Minutes into run')
ylabel('Doppler residuals [m]')

end



























