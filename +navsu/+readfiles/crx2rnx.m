function [status, sysout] = crx2rnx(filename)
% crx2rnx
% DESCRIPTION:
%   Decompress Hatanaka compressed RINEX observation files (.**d, .crx) 
%   using crx2rnx.exe 
% INPUT:
%   filename   - file to be decompressed
% OUTPUT:
%   status     - status output from crx2rnx.exe
%
% See also: navsu.readfiles.loadRinexObs 


% get file path
s = what('+navsu');
thirdpartypath = [s.path '\+thirdparty\'];

crx2rnxLoc = [thirdpartypath '\crx2rnx.exe'];
% Open desired file using windows
[status, sysout] = system(['"' crx2rnxLoc '" "' filename '" -f']);
end