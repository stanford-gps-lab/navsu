function [status,result] = unzipFile(filename,outLocation)
% unzipFile
% DESCRIPTION:
%   Unzip a given file using the portable version of 7zip included in the
%   repo if this is a windows machine
% INPUT:
%   filename    - name of the file to be parsed
% OPTIONAL INPUTS:
%   outLocation - If desired, you can specify a different output location
%                 than the same folder as the input file
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

if ispc
    % Windows machines use the included copy of 7zip
    s = what('+navsu');
    
    utilityPath = fullfile(s.path, '/+thirdparty/7zip/');
    
    loc7zip = fullfile(utilityPath ,'7za.exe');
    % Use 7zip to open!
    [status,result] = system(['"' loc7zip '" -y x ' '"' filename '"' ' -o' '"' outLocation '\"']);
    
else
    
    if endsWith(filename, '.zip', 'IgnoreCase', true), % {gunzip, uncompress}: these fail (tested on macOS 10.14)
        
        try
            result = unzip(filename, outLocation); status = 0;            
        catch
            result = []; status = 1;
        end
        
    elseif endsWith(filename, '.gz', 'IgnoreCase', true),  % {unzip, uncompress}: these fail (tested on macOS 10.14)
        
        % MATLAB built-in gunzip introduced before R2006a)
        try
            result = gunzip(filename, outLocation); status = 0;
        catch
            result = []; status = 1;
        end
        %   Alternatively, try host machine's native gunzip executable:
        %   [status,result] = system(['gunzip ' filename]);
        
                
    elseif endsWith(filename, '.Z'),
        
        % host machine's native executable
        [status,result] = system(['gunzip ' filename]);
        
        % (MAY ALSO WORK:) host machine's native uncompress executable, if available
        % [status,result] = system(['uncompress ' filename]);
        
        % (DOES NOT WORK:) as of R2020a, MATLAB's built-in gunzip version is older
        % than that bundled with recent versions of macOS (tested on 10.14)
        % unzippedFilename = gunzip(filename, outLocation);
        
    else
        
        error('Sorry, nothing is currently implemented to unzip files for non-windows machines. See <https://www.7-zip.org/download.html> for possible alternatives.');
        
    end
    
end


end