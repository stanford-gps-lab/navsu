function prn = svn2prn(svn, epoch, const)
% prn2svn
% DESCRIPTION:
%   Map from SVN to PRN at a given time.
%
% INPUT:
%   svn     - satellite vehicle number
%   epoch   - time is in fractional year (e.g 2012.233 very crudely) or 
%             Julian day (e.g. 2444239.5)
%   
% OUTPUT:
%   prn     - PRN
%
% See also: navsu.svprn.prn2svn, navsu.time.cal2jd,
%           navsu.svprn.constSvnData

% PRN is column 2, SVN is column 1:
prn = navsu.svprn.mapSatData(2, 1, svn, epoch, const);

end