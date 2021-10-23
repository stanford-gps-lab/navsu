function files = curlGetDirectoryContents(site, netrcFile, cookieFile)
% Pull directory contents using curl

% make sure the path is ending with a slash
if ~strcmp(site(end), '/')
    site = [site, '/'];
end

if ispc && nargin == 3 && isfile(netrcFile) && isfile(cookieFile)
    curlCall = ['curl --silent -c "' cookieFile ...
                '" -n --netrc-file "' netrcFile '" -L "' ...
                site '*?list"'];
else
    curlCall = ['curl -u anonymous:fabianr@stanford.edu ' ...
                '--ftp-ssl ' site];
end

[curlFeedbackCode, output] = system(curlCall);

% gracefully warn user of failure
if curlFeedbackCode > 0
    warning(['Failed to access %s\n', ...
             'Could not retrieve ftp directory listing.\n', ...
             'cURL exited with code %i.'], site, curlFeedbackCode);
end

if ispc
    % legacy code. Has to be tested on windows!!
    files = textscan(output, '%s%f');
    files = files{1};
else
    % this works with the output of the ftp call
    byLineOutput = split(output, newline);
    
    useLines = find(~cellfun(@isempty, byLineOutput));
    nLines = length(useLines);
    
    files = cell(nLines, 1);
    
    for lineId = 1:nLines
        s = textscan(byLineOutput{useLines(lineId)}, '%s');
        files{lineId} = s{1}{end};
    end
end


end