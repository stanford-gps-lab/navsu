function [status,result] = unzipFile(filename,outLocation)

% No output location specified- unzip to same directory
if nargin == 1 
   outLocation = fileparts(filename); 
end
s = what('+utility');

utilityPath = [s.path '/+thirdparty/7zip/'];

loc7zip = [utilityPath '7za.exe'];
    % Use 7zip to open!
    [status,result] = system(['"' loc7zip '" -y x ' '"' filename '"' ' -o' '"' outLocation '"']);
% end

end