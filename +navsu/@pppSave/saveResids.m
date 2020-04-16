function saveResids(obj,measMat,resids,epoch,el,prnConstInds)

% Save the range residuals
[~,indsSave] = ismember(measMat(measMat(:,end) == 1 | ...
    measMat(:,end) == 2,[1 2 3 6]),[obj.gnssData.range.PRN(:) obj.gnssData.range.constInds(:) ...
    obj.gnssData.range.sig(:) obj.gnssData.range.ind(:)],'rows');
residsSavei = nan(size(obj.gnssData.range.resids,1),size(obj.gnssData.range.resids,2));
residsSavei(indsSave) = resids(measMat(:,6) == 1 | measMat(:,6) == 2);


indEpoch = find(obj.gnssData.epochs == epoch);
obj.gnssData.range.resids(:,:,indEpoch) = residsSavei;

% Save the doppler residuals
[~,indsSave] = ismember(measMat(measMat(:,end) == 3 ,[1 2 3]),...
    [obj.gnssData.doppler.PRN(:) obj.gnssData.doppler.constInds(:) ...
    obj.gnssData.doppler.sig(:) ],'rows');
residsSavei = nan(size(obj.gnssData.doppler.resids,1),size(obj.gnssData.doppler.resids,2));
residsSavei(indsSave) = resids(measMat(:,6) == 3);

indEpoch = find(obj.gnssData.epochs == epoch);
obj.gnssData.doppler.resids(:,:,indEpoch) = residsSavei;

% Save elevation per satellite
[~,indsEl] = ismember(prnConstInds,[obj.gnssData.PRN' obj.gnssData.constInds'],'rows');
obj.gnssData.el(indsEl,indEpoch) = el;


end