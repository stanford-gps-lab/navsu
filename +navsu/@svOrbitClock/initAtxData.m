function initAtxData(obj,filenameAtx)

atxData = utility.readfiles.readAtx(filenameAtx);

% remove data that isn't for satellites
atxData = atxData([atxData(:).type] ~= 0);

obj.atx = atxData;

end