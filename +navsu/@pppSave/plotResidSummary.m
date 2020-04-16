function plotResidSummary(obj,varargin)
%% Randon optional inputs
p = inputParser;

p.addParameter('dcbCorr',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
dcbCorr   = res.dcbCorr;  

%% compute mean and std for all measurements

means = nanmean(obj.gnssData.range.resids,3);
stds = nanstd(obj.gnssData.range.resids,0,3);

indsPr = find(obj.gnssData.range.ind(:,1) == 1 & any(~isnan(means)')');
indsPh = find(obj.gnssData.range.ind(:,1) == 2 & any(~isnan(means)')');

prns = obj.gnssData.PRN(1,:)';

% first populated PRN
indsAvail = find(~cellfun(@isempty,obj.gnssData.range.rnxCode(1,:)));

%% Plot mean and standard deviations for each PRN
figi = obj.figResidsSummary;
if ~isvalid(figi)
    figi = figure;
end
figi.Visible = 'on';
% figi.Position = [figi.Position(1) figi.Position(2)+(figi.Position(4)-688) 932         688];
clf;
figi.Position = [figi.Position(1:2)-630/2 800 630];
ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.09 0.05]);

colors = lines(length(indsPr));

axes(ha(1));
for idx = 1:length(indsPr)
    errorbar(1:size(means,2),means(indsPr(idx),:),stds(indsPr(idx),:),'o','color',colors(idx,:))
    
    hold on;
    
    if ~isempty(dcbCorr)
        plot(dcbCorr(indsPr(idx),:),'^','color',colors(idx,:),'markersize',5)
    end
end
grid on
% ylim([-10 10])
ylabel('Code residuals')
% legend(obj.gnssData.rnxCode(indsPr,indsAvail(2)))

xticks(1:length(prns))
xticklabs = repmat('  ',length(prns),1);
xlim([0.5 length(prns)+0.5])
xticklabels(xticklabs)
title('Mean and std of code and carrier residuals')

axes(ha(2));
for idx = 1:length(indsPh)
    errorbar(1:size(means,2),means(indsPh(idx),:),stds(indsPh(idx),:),'o')
    hold on;
end
grid on
ylim([-0.1 0.1])
ylabel('Carrier residuals')

xticks(1:length(prns))
xticklabs = repmat('  ',length(prns),1);
xticklabels(xticklabs)
xlim([0.5 length(prns)+0.5])

% Plot doppler residuals

meansDop = nanmean(obj.gnssData.doppler.resids,3);
stdsDop = nanstd(obj.gnssData.doppler.resids,0,3);
axes(ha(3));
for idx = 1:size(meansDop,1)
    errorbar(1:size(meansDop,2),meansDop(idx,:),stdsDop(idx,:),'o')
    hold on;
end
grid on
xlabel('PRN')
ylabel('Doppler residuals')
% legend(obj.gnssData.rnxCode(indsPh,indsAvail(2)))

xticks(1:length(prns))
xticklabs = repmat('  ',length(prns),1);
xticklabs(1:2:end,:) = [num2str(prns(1:2:end))];
xticklabels(xticklabs)
xlim([0.5 length(prns)+0.5])

%% Plot GLONASS residuals vs frequency offset
if 0
if any(any(obj.gnssData.constInds == 2))
    figi = obj.figResidsSummaryGlo;
    if ~isvalid(figi)
        figi = figure;
    end
    figi.Visible = 'on';
    freq0 = 1602000000;
    for idx = 1:length(indsPr)
        %     if idx == 2
        %         errorbar(1:size(means,2),means(indsPr(idx),:),stds(indsPr(idx),:),'o','color',colors(idx,:))
        %     else
        indsGlo = find(obj.gnssData.range.constInds(indsPr(idx),:) == 2);
        errorbar(obj.gnssData.range.freqs(indsPr(idx),indsGlo)-freq0,means(indsPr(idx),indsGlo),stds(indsPr(idx),indsGlo),'o','color',colors(idx,:))
        %     end
        hold on;
        
        if ~isempty(dcbCorr)
            plot(dcbCorr(indsPr(idx),:),'^','color',colors(idx,:),'markersize',5)
        end
    end
    
end

end


end

























