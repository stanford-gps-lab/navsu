function curlDownloadDirectory(file,localDir,netrcFile,cookieFile)
% Use curl to download a directory

% pull the filename
[remoteDir, ~, ~] = fileparts(file);

% get list of all files
fileNames = navsu.ftp.curlGetDirectoryContents(remoteDir, netrcFile, cookieFile);
filePaths = fullfile(remoteDir, fileNames);


if ispc
    % Functionality to download all at once using curl- much faster than
    % downloading individual files. 
    
    % check if the local folder exists
    if ~exist(localDir,'dir')
        mkdir(localDir);
        % download each file
    end
    
    filenameAny = ['*'];
    
    [~,output] = system(['cd "' localDir '" & curl -c "' cookieFile '" --silent -n --netrc-file "' netrcFile '" -L -O --remote-name-all "' remoteDir '/' filenameAny '" ']);
    
    % unzip the downloaded files
    navsu.readfiles.unzipFile(fullfile([localDir '_']),localDir,true);
    
    % delete the original compressed file
    delete(fullfile([localDir '_']));
    
else
    % download each file
    for fI = 1:length(filePaths)
        
        navsu.ftp.curlDownloadSingleFile(filePaths{fI},localDir, netrcFile, cookieFile);
        
    end
    
    % unzip the downloaded files
    navsu.readfiles.unzipFile(fullfile([localDir '_']))
    
    % delete the original compressed file
    
    delete(fullfile([localDir '_']))
    
end

end