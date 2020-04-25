function fullMeas = syncMeas(varargin)

measUngrouped = cat(1,varargin{:});

nMeasUngrouped = length(measUngrouped);

epochsUngrouped = nan(nMeasUngrouped,1);

for idx = 1:nMeasUngrouped
   epochsUngrouped(idx) = measUngrouped{idx}.epochs; 
end

epochsUnique = unique(epochsUngrouped);

fullMeas = cell(length(epochsUnique),1);

for idx =1:length(epochsUnique)
    fullMeas{idx} = measUngrouped(epochsUngrouped == epochsUnique(idx));
end


end