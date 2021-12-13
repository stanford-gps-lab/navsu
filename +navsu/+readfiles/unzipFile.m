function [status,result] = unzipFile(filename, outLocation, overwrite)
% unzipFile
% DESCRIPTION:
%   Unzip a given file using the portable version of 7zip included in the
%   repo if this is a windows machine. Can alternatively be called on a
%   folder to unpack all contained archives. Recurses on subfolders.
% INPUT:
%   filename    - name of the file or folder to be parsed
% OPTIONAL INPUTS:
%   outLocation - If desired, you can specify a different output location
%                 than the same folder as the input file
%   overwrite   - Should local files be overwritten?
%
% OUTPUT:
%   status      - status output by 7zip or other uncompression utility
%   result      - output of uncompression command, as sent to stdout
% 
% See also: 

% No output location specified- unzip to same directory
if nargin == 1 
   outLocation = fileparts(filename); 
end
if nargin < 3
    overwrite = false;
end

% check if file already exists locally
[~, fileName, ~] = fileparts(filename);
if isfile(fullfile(outLocation, fileName)) && ~overwrite
    % unzipped file already exists!
    result = []; status = 0;
    return
end

if ispc
    % Windows machines use the included copy of 7zip
    s = what('+navsu');
    
    utilityPath = fullfile(s.path, '/+thirdparty/7zip/');
    
    loc7zip = fullfile(utilityPath ,'7za.exe');
    % Use 7zip to open!
    [status,result] = system(['"' loc7zip '" -y x ' '"' filename '"' ' -o' '"' outLocation '\"']);
    
else
    
    if endsWith(filename, '.zip', 'IgnoreCase', true) % {gunzip, uncompress}: these fail (tested on macOS 10.14)
        
        try
            result = unzip(filename, outLocation); status = 0;            
        catch
            result = []; status = 1;
        end
        
    elseif endsWith(filename, '.gz', 'IgnoreCase', true)  % {unzip, uncompress}: these fail (tested on macOS 10.14)
        
        % MATLAB built-in gunzip introduced before R2006a)
        try
            result = gunzip(filename, outLocation); status = 0;
        catch
            result = []; status = 1;
        end
        %   Alternatively, try host machine's native gunzip executable:
        %   [status,result] = system(['gunzip ' filename]);
        
                
    elseif endsWith(filename, '.Z')
        
        % host machine's native executable
        [status, result] = system(['gunzip "' filename '"']);
        
        % (MAY ALSO WORK:) host machine's native uncompress executable, if available
        % [status,result] = system(['uncompress ' filename]);
        
        % (DOES NOT WORK:) as of R2020a, MATLAB's built-in gunzip version is older
        % than that bundled with recent versions of macOS (tested on 10.14)
        % unzippedFilename = gunzip(filename, outLocation);
        
    elseif isfolder(filename)
        % recursively unzip everything in the folder
        d = dir(filename);
        
        status = ones(length(d), 1); % 0 = success, 1 = failure
        result = cell(length(d), 1);
        % to avoid stack overflows from "." and ".."
        toUnpack = ~arrayfun(@(x) startsWith(x.name, '.') ...
                                || endsWith(x.name, '.'), d);
        for dI = find(toUnpack')
            % try unzipping every file
            [status(dI), result{dI}] = navsu.readfiles.unzipFile( ...
                fullfile(filename, d(dI).name), outLocation, overwrite);
        end
        % compile output variables
        status =  max(status);
        haveResult = ~cellfun(@isempty, result);
        if any(haveResult)
            result = join([result{haveResult}], '\n');
        else
            result = [];
        end
    else
        result = []; status = 1;
        fprintf('Missing the means to unpack %s on non-windows machines.\n', filename);
%         error(['Sorry, nothing is currently implemented to unzip these files for non-windows machines. ', ...
%                'See <https://www.7-zip.org/download.html> for possible alternatives.']);
        
    end
    
end


end