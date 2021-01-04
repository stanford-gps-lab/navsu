function curlDownloadSingleFile(file,localDir,netrcFile,cookieFile)
% Use curl to download a single file

% pull the filename
[remoteDir,filename,ext] = fileparts(file);

filename = [filename ext];

% check if the local folder exists
if ~exist(localDir,'dir')
    mkdir(localDir);
end

[~,output] = system(['curl -c "' cookieFile '" --silent -n --netrc-file "' netrcFile ' " -L -o "' fullfile(localDir, filename) '" "' remoteDir '/' filename '" ']);


end