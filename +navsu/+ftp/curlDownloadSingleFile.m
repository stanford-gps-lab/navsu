function curlDownloadSingleFile(remoteFile, localDir, varargin)
% Use curl to download a single file from a remote directory.
% 
%   curlDownloadSingleFile(remoteFile, localDir, netrcFile, cookieFile)
%   curlDownloadSingleFile(remoteFile, localDir)
%   
%   Can be called with and without netrc and cookie file. On windows opts
%   for faster http download if netrc and cookie file are supplied. See:
%   https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
%   https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html


if startsWith(localDir, '~')
    warning('Local directory file path cannot start with "~" but must be fully specified!')
end

% pull the filename
[~, filename, ext] = fileparts(remoteFile);

localFile = fullfile(localDir, [filename ext]);

% check if the local folder exists
if ~exist(localDir, 'dir')
    mkdir(localDir);
end

if navsu.ftp.curlUseHttp(varargin{:})
    % use http protocol
    system(['curl --silent -c "' varargin{2} ...
            '" -n --netrc-file "' varargin{1} ...
            ' " -L -o "' localFile '" "' remoteFile '" ']);
elseif ispc
    % can use ftp-ssl protocol, but need to force timeout
    system(['curl --silent --speed-time 1 --speed-limit 10 ' ...
            '-u anonymous:fabianr@stanford.edu ' ...
            '-o "' localFile '" --ftp-ssl "' remoteFile '"']);
else
    % use ftp-ssl protocol
    system(['curl --silent ' ...
            '-u anonymous:fabianr@stanford.edu ' ...
            '-o "' localFile '" --ftp-ssl "' remoteFile '"']);
end

% fail gracefully: warn user of failed download
if ~isfile(localFile)
    warning(['cURL download of ', remoteFile, ' failed!']);
end

end