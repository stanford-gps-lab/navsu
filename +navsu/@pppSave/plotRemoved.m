function plotRemoved(obj)


figi = obj.figRemoved;
if ~isvalid(figi)
    figi = figure;
end
figi.Visible = 'on';
clf; hold on;
%%
measRemoved = obj.measRemoved(~isnan(obj.epochsRemoved),:);
epochsRemoved = obj.epochsRemoved(~isnan(obj.epochsRemoved));

prns = obj.gnssData.PRN(1,:)';
constInds = obj.gnssData.constInds(1,:)';
sInds = (1:length(prns))';
epochs = obj.gnssData.epochs;

[~,yInds] = ismember(measRemoved(:,1:2),[prns constInds],'rows');
[~,xInds] =  ismember(epochsRemoved,epochs);

yPlotMat = repmat(sInds,1,length(epochs));
xPlotMat = repmat(1:length(epochs),length(sInds),1);

markers = {'g.','c.','ro','r.','b^','k'};
markerSize = [5 2 5 8 2 2];
legText = {'Meas Used','Low Elevation','Code Removed','Carr Removed','Cycle Slip','Number of satellites used'};
for idx = 1:length(legText)
    plot(-10,-10,markers{idx},'markerSize',markerSize(idx))
end

for idx = 1:length(legText)
    switch idx
        case 1 % Meas used
            linInd = find(squeeze(sum(~isnan(obj.gnssData.range.resids),1)) > 0);
            
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