function initBroadcastData(obj,years,doys,varargin)


p = inputParser;
p.addParameter('DOWNLOAD',true);   

% parse the results
parse(p, varargin{:});
res = p.Results;
DOWNLOAD        = res.DOWNLOAD;        % indicator to check for downloads and download

if numel(years) == 1
    % bring it to the same size as doys
    years = years * ones(size(doys));
end
if numel(years) ~= numel(doys)
    error('years and doys inputs must have same number of elements.');
end

if DOWNLOAD   
    navsu.ftp.download(16, years, doys, obj.settings);
end

% Load the data!

% if sum(obj.settings.constUse) > 1
    % multi-const
    fileEnding = 'p';
% else
%     % GPS only
%     fileEnding = 'n';
% end

for dayi = 1:length(doys)
    yr = years(dayi);
    doy = doys(dayi);
    % set file path, find broadcast file within
    filePath = fullfile(obj.settings.navMgxDir, ...
                        num2str(yr), ...
                        num2str(doy, '%03d'));
    folderDir = dir(filePath);
    files = arrayfun(@(x) x.name, folderDir, 'UniformOutput', false);
    fileId = cellfun(@(fName) endsWith(lower(fName), ...
                     ['.', num2str(mod(yr, 100)), fileEnding]), files);
    if ~any(fileId)
        % try MGEX file
        fileId = cellfun(@(fName) endsWith(fName, '.rnx'), files);
    end
    if ~any(fileId)
        % try GPS only file
        fileId = cellfun(@(fName) endsWith(lower(fName), ...
                     ['.', num2str(mod(yr, 100)), 'n']), files);
    end
    fileName = files{fileId};
        
    % now parse the file
    constCell = num2cell(obj.settings.constUse);
    ephi = navsu.readfiles.loadRinexNav(fullfile(filePath, fileName), ...
        'constellations', navsu.readfiles.initConstellation(constCell{:}));

    % concatenate structs from multiple days
    if dayi == 1
        eph = ephi;
    else
        eph = mergeEphStructsRecursive(eph, ephi);
    end
    
end

% Put it in the object
obj.BEph = eph;


% internal helper function to merge the structs:

function s = mergeEphStructsRecursive(s1, s2)
    % Merge two eph structures. Vertically concatenates all elements.
    % 
    % s = mergeEphStructsRecursive(obj, s1, s2)
    % 
    % Recursively works on fields that are structures themselves.
    % Both structures need to have the same architecture.
    % 

    % initialize as equal to s1
    s = s1;

    % now append all fields from s2, respect arrays
    nS = length(s2);
    for iS2 = 1:nS

        f2 = fields(s2(iS2));
        for if2 = 1:length(f2)
            if isstruct(s2(iS2).(f2{if2}))
                % recursive case
                s(iS2).(f2{if2}) = mergeEphStructsRecursive(s1(iS2).(f2{if2}), s2(iS2).(f2{if2}));
            elseif isfield(s1(iS2), f2{if2})
                % base case 1: s1 has the field as well
                s(iS2).(f2{if2}) = vertcat(s1(iS2).(f2{if2}), s2(iS2).(f2{if2}));
            else
                % base case 2: s1 does not yet have the field
                s(iS2).(f2{if2}) = s2(iS2).(f2{if2});
            end

        end
    end
end

end

