function svn = prn2svn(prn, epoch, varargin)
% prn2svn(prn, epoch, const, source)
% DESCRIPTION:
%   Map from PRN to SVN at a given time.
%
% INPUT:
%   prn     - N x 1 vector of PRNs
%   epoch   - time is in fractional year (e.g 2012.233 very crudely) or 
%             Julian day (e.g. 2444239.5)
%   const   - OPTIONAL N x 1 vector of constellation indices (defaults to
%             GPS)
%   source  - OPTIONAL data source (only applies to GLONASS)
%   
% OUTPUT:
%   svn     - satellite vehicle number.
%
% See also: navsu.svprn.prn2x, navsu.svprn.svn2prn, navsu.time.cal2jd

svn = navsu.svprn.prn2x(1, prn, epoch, varargin{:});

end