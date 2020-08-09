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


% get file path -- the use of 'filesep' isn't strictly necessary here since
% this is implicitly directed to Windows machines, but should cause no issues
% and may be more future-proof (e.g. if host is a *nix machine, but runs Wine
% or similar to handle .exe files).
s = what('+navsu');
thirdpartypath = [s.path filesep '+thirdparty' filesep];
crx2rnxLoc = [thirdpartypath 'crx2rnx.exe'];

% Open desired file using windows
[status, sysout] = system(['"' crx2rnxLoc '" "' filename '" -f']);
end