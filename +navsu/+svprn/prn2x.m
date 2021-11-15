function outVal = prn2x(colId, prn, epoch, varargin)
% prn2x(colId, prn, epoch, const, source)
% DESCRIPTION:
%   Map from PRN to satellite specific value at a given time.
%
% INPUT:
%   colId   - index of column of desired data from the svn table (see
%             navsu.svprn.constSvnData)
%   prn     - N x 1 vector of queried PRN numbers
%   epoch   - time is in fractional year (e.g 2012.233 very crudely) or 
%             Julian day (e.g. 2444239.5)
%   const   - OPTIONAL constellation indices, integer or vector of same
%             size as prn input (defaults to GPS)
%   source  - OPTIONAL data source (only applies to GLONASS)
%   
% OUTPUT:
%   outVal  - desired output value of same dimensions as prn input
%
% See also: navsu.svprn.prn2svn, navsu.svprn.prn2FreqChanGlonass

% PRNs are stored in column 2
outVal = navsu.svprn.mapSatData(colId, 2, prn, epoch, varargin{:});

end