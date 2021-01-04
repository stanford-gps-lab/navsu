function curlDownloadDirectory(file,localDir,netrcFile,cookieFile)
% Use curl to download a single file

% pull the filename
[remoteDir,filename,ext] = fileparts(file);

filename = [filename ext];

% check if the local folder exists
if ~exist(localDir,'dir')
    mkdir(localDir);
end



filenameAny = ['*'];

[~,output] = system(['cd "' localDir '" & curl -c "' cookieFile '" --silent -n --netrc-file "' netrcFile '" -L -O --remote-name-all "' remoteDir '/' filenameAny '" ']);


% unzip the downloaded file
navsu.readfiles.unzipFile(fullfile([localDir '_']))

% delete the original compressed file

delete(fullfile([localDir '_']))


end