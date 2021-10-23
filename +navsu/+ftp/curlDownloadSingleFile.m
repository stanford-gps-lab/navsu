function curlDownloadSingleFile(file, localDir, netrcFile, cookieFile)
% Use curl to download a single file

if startsWith(localDir, '~')
    warning('Local directory file path cannot start with "~" but must be fully specified!')
end

% pull the filename
[remoteDir, filename, ext] = fileparts(file);

filename = [filename ext];

localFile = fullfile(localDir, filename);

% check if the local folder exists
if ~exist(localDir,'dir')
    mkdir(localDir);
end

if ispc && nargin == 4 && isfile(netrcFile) && isfile(cookieFile)
    system(['curl --silent -c "' cookieFile ...
            '" -n --netrc-file "' netrcFile ...
            ' " -L -o "' localFile ...
            '" "' remoteDir '/' filename '" ']);
else
    % use ftp-ssl protocol on mac:
    system(['curl --silent -u anonymous:fabianr@stanford.edu -o "' ...
            localFile '" --ftp-ssl "' remoteDir '/' filename '"']);
end

% fail gracefully: warn user of failed download
if ~isfile(localFile)
    warning(['cURL download of ', filename, ' from ', remoteDir, ' failed!']);
end

end