function initAtxData(obj,filenameAtx)
% Download (if necessary) and load IGS antenna phase center (.atx) file
% DESCRIPTION:
%   Initialize the IGS antenna phase center file
%   
% INPUT:
%   filenameAtx - name and path of the IGS antenna phase center file (.atx)
%  
% OUTPUT:
%   The object will have a brand new antenna phase center structure!  
%
% See also: 

atxData = navsu.readfiles.readAtx(filenameAtx);

% remove data that isn't for satellites
atxData = atxData([atxData(:).type] ~= 0);

obj.atx = atxData;

end