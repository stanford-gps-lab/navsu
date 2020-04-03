function initAtxData(obj,filenameAtx)

atxData = ReadATX(filenameAtx);
% Find antenna data for the receiver of interest
% atxRx = atxData( find(~cellfun(@isempty,strfind({atxData.block},antMod))));

% remove data that isn't for satellites
atxData = atxData([atxData(:).type] ~= 0);

obj.atx = atxData;

end