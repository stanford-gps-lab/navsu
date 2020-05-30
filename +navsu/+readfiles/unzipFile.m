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
%   status      - status output by 7zip
%
% See also: 

% No output location specified- unzip to same directory
if nargin == 1 
   outLocation = fileparts(filename); 
end

if ispc
    % Windows machines use the included copy of 7zip
    s = what('+navsu');
    
    utilityPath = [s.path '/+thirdparty/7zip/'];
    
    loc7zip = [utilityPath '7za.exe'];
    % Use 7zip to open!
    [status,result] = system(['"' loc7zip '" -y x ' '"' filename '"' ' -o' '"' outLocation '"']);
    
else
    error('Sorry, nothing is currently implemented to unzip files for non-windows machines')
    
    % TODO:
    % Something should be implemented here to uncompress files for
    % mac/linux machines.  7zip happily unzips many different compression
    % types (.gz, .zip, .Z) , so whatever is done here should do the same.  
    
        
end


end