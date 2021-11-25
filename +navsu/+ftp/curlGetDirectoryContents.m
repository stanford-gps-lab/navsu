function files = curlGetDirectoryContents(remoteDir, varargin)
% Pull directory contents using curl. Returns a cell array of strings. Each
% string specifies one file in the remote directory.
%    
%   files = curlGetDirectoryContents(site, netrcFile, cookieFile)
%   files = curlGetDirectoryContents(site)
%   
%   Can be called with netrc and cookie file on windows for slightly faster
%   response using http instead of ftp. See:
%   https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
%   https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html



% make sure the path is ending with a slash
if ~endsWith(remoteDir, '/')
    remoteDir = [remoteDir, '/'];
end

% use ftp or http?
useHttp = navsu.ftp.curlUseHttp(varargin{:});

if useHttp
    curlCall = ['curl --silent -c "' varargin{2} ...
                '" -n --netrc-file "' varargin{1} ...
                '" -L "' remoteDir '*?list"'];
elseif ispc
    % need to call ftp protocol with a timeout
    curlCall = ['curl --silent --speed-time 1 --speed-limit 10 '...
                '-u anonymous:fabianr@stanford.edu ' ...
                '--ftp-ssl ' remoteDir];
else
    % use ftp protocol on the mac
    curlCall = ['curl --silent ' ...
                '-u anonymous:fabianr@stanford.edu ' ...
                '--ftp-ssl ' remoteDir];
end

[curlFeedbackCode, output] = system(curlCall);

% gracefully warn user of failure (unless normal timeout code for ftp on pc)
if curlFeedbackCode > 0 && (useHttp || ~ispc)
    warning(['Failed to access %s\n', ...
             'Could not retrieve ftp directory listing.\n', ...
             'cURL exited with code %i.'], remoteDir, curlFeedbackCode);
end

% now parse the output of the curl command
if useHttp
    % this works with the output of the http call
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