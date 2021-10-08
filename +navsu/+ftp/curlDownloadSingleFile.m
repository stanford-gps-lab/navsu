function curlDownloadSingleFile(file, localDir, netrcFile, cookieFile)
% Use curl to download a single file

% pull the filename
[remoteDir, filename, ext] = fileparts(file);

filename = [filename ext];

% check if the local folder exists
if ~exist(localDir,'dir')
    mkdir(localDir);
end

if ispc
    system(['curl --silent -c "' cookieFile ...
            '" -n --netrc-file "' netrcFile ...
            ' " -L -o "' fullfile(localDir, filename) ...
            '" "' remoteDir '/' filename '" ']);
else
    % use ftp-ssl protocol on mac:
    system(['curl --silent -u anonymous:fabianr@stanford.edu -o "' ...
            fullfile(localDir, filename) ...
            '" --ftp-ssl "' remoteDir '/' filename '"']);
end

end