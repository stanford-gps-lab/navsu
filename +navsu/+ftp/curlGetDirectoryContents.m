function files = curlGetDirectoryContents(site,netrcFile,cookieFile)
% Pull directory contents using curl

[~,output] = system([' curl -c ' cookieFile ' --silent -n --netrc-file "' netrcFile '" -L "' ...
     site '*?list"']);

files = textscan(output,'%s%f');

files = files{1};


end