function indsFind = strFindCell(str,pattern,noFind,exact)
% strFindCell
% DESCRIPTION:
%   Searches for a string patten in a cell array of strings
% INPUT:
%   str     - cell array of strings
%   pattern - string to search for
% OPTIONAL INPUTS:
%   noFind  - Returns true at each index when string was found, otherwise
%             return the index numbers where the string was found
%   exact   - Strings in the cell array must match exactly.  Otherwise,
%             searches for substrings. 
% OUTPUT:
%   indsFind - index numbers where the string was found in the cell array
%
% See also: 

if nargin < 3
    noFind = 0;
end

if nargin < 4
    exact = 0;
end

if ~exact
    if noFind
        indsFind = ~cellfun(@isempty, strfind(str,pattern));
    else
        indsFind = find(~cellfun(@isempty, strfind(str,pattern)));
    end
else
    if noFind
        indsFind = strcmp(str,pattern);
    else
        indsFind = find(strcmp(str,pattern));
    end
end


end