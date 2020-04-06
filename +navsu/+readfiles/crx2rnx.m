function [status, sysout] = crx2rnx(filename)
% Wrapper function for using crx2rnx (Hatanaka de-compression .exe)
% [pathstr,name,ext] = fileparts(filename);
% need to change name of file to .crx
% if ~strcmp(ext,'.crx')
%     copyfile(filename,[pathstr '\' name '.crx'])
% end

% get file path
s = what('+navsu');
thirdpartypath = [s.path '\+thirdparty\'];

crx2rnxLoc = [thirdpartypath '\crx2rnx.exe'];
% Open desired file using windows
[status, sysout] = system(['"' crx2rnxLoc '" "' filename '" -f']);
end