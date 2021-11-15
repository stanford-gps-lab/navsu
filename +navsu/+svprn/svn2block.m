function [block, blockNums] = svn2block(svn,const,flag_number,source)
% Function to translate GNSS SVN to block 
% CONST is index of GNSS constellation:
% GPS = 1, GLO = 2, GAL = 3, BDS = 4
%   If not included, defaults to GPS

if nargin < 2
    const = 1;
end
% Check for string constellation type input and convert to number
if ischar(const)
    constNames = {'GPS','GLO','GAL','BDS'};
    const = find(~cellfun(@isempty,(strfind(constNames,const))));
    if isempty(const)
        warning('Invalid constellation name- defaulting to GPS')
        const = 1;
    end
end

if nargin < 3
    flag_number = 0;
end

% GLONASS database source
if nargin < 4
    if const == 2
        source = 3;
    else
        source = 1;
    end
end

[svndata, blockText] = navsu.svprn.constSvnData(const,source);

% preallocate outputs
nSvn = length(svn);
blockNums = NaN(1, nSvn);
block = repmat({''}, 1, nSvn);

for i = 1:nSvn
    idx = svndata(svndata(:,1) == svn(i), 13);
    %     if flag_number
    if ~isempty(idx)
        blockNums(i) = idx(1);
    end
    %     else
    if ~isempty(idx)
        block{i} = blockText{idx};
    end
    %     end
end

if flag_number
    block = blockNums;
end
end











