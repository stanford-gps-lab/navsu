function [obsMap,epochs] = syncMeas(obs)

epochsFull =[];
for idx = 1:length(obs)
    epochsFull = [epochsFull; obs{idx}.epochs];
end
epochs = unique(epochsFull);
obsMap = nan(length(epochs),length(obs));
for idx = 1:length(obs)
    [~,ixb] = ismember(epochs,obs{idx}.epochs);
    obsMap(:,idx) = ixb;
end

end