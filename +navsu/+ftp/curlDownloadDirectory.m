function curlDownloadDirectory(file, localDir, netrcFile, cookieFile)
% Use curl to download a directory

% pull the filename
[remoteDir, ~, ~] = fileparts(file);

if startsWith(localDir, '~')
    warning('Local directory file path cannot start with "~" but must be fully specified!')
end

if ispc && nargin == 4 && isfile(netrcFile) && isfile(cookieFile)
    % Functionality to download all at once using curl- much faster than
    % downloading individual files. But requires netrc and cookie file.
    % This currently only works on windows.
    
    % check if the local folder exists
    if ~exist(localDir, 'dir')
        mkdir(localDir);
    end
    
    % download each file
    [~,output] = system(['cd "' localDir '" & curl -c "' cookieFile ...
                         '" --silent -n --netrc-file "' netrcFile ...
                         '" -L -O --remote-name-all "' remoteDir '/*" ']);
    
    % unzip the downloaded files
    navsu.readfiles.unzipFile(fullfile([localDir '_']), localDir, true);
    
    % delete the original compressed file
    delete(fullfile([localDir '_']));
    
else
    % loop over all files, pull them one by one
    
    % get list of all files
    fileNames = navsu.ftp.curlGetDirectoryContents(remoteDir, netrcFile, cookieFile);
    filePaths = fullfile(remoteDir, fileNames);

    % download and unzip each file, keep only uncompressed version
    for fI = 1:length(filePaths)

        navsu.ftp.curlDownloadSingleFile(filePaths{fI}, localDir, netrcFile, cookieFile);

        % unzip the downloaded file
        navsu.readfiles.unzipFile(fullfile(localDir, fileNames{fI}));
        % delete the original compressed file
        delete(fullfile(localDir, fileNames{fI}));
    end

end

end