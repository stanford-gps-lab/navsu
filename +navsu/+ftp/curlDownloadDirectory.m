function curlDownloadDirectory(file,localDir,netrcFile,cookieFile)
% Use curl to download a directory

% pull the filename
[remoteDir, ~, ~] = fileparts(file);

% get list of all files
fileNames = navsu.ftp.curlGetDirectoryContents(remoteDir, netrcFile, cookieFile);
filePaths = fullfile(remoteDir, fileNames);

% download each file
for fI = 1:length(filePaths)
    
    navsu.ftp.curlDownloadSingleFile(filePaths{fI}, netrcFile, cookieFile);
    
end

% unzip the downloaded files
navsu.readfiles.unzipFile(fullfile([localDir '_']))

% delete the original compressed file

delete(fullfile([localDir '_']))

end