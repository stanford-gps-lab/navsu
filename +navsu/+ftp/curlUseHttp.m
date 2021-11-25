function useHttp = curlUseHttp(varargin)
%Determines whether or not the http or ftp axess to the cddis server should
%be used. Opts for http on windows, if a valid cookie and netrc file are
%passed as inputs. Returns false unless it is called with two character
%vectors as inputs that both represent valid file paths.
%   
%   useHttp = navsu.ftp.curlUseHttp(netrcFile, cookieFile)

useHttp = ispc && numel(varargin) == 2 ...
               && all(cellfun(@ischar, varargin)) ...
               && all(cellfun(@isfile, varargin));
end