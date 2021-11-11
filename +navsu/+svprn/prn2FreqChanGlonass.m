function freqChan = prn2FreqChanGlonass(prn, epoch, source)
% prn2FreqChanGlonass
% DESCRIPTION:
%   Map from GLONASS PRN to the FDMA frequency channel at a given time.
% INPUT:
%   prn     - PRN
%   epoch   - time is in fractional year (e.g 2012.233 very crudely) or 
%             Julian day (e.g. 2444239.5)
%   
% OUTPUT:
%   freqChan - GLONASS FDMA frequency channel (-7 to 7) 
%
% See also: navsu.svprn.prn2svn, navsu.svprn.svn2prn, navsu.time.cal2jd

if nargin < 2
    error('You must specify an svn and a time');
end

if nargin < 3
    % Source of database (only applies to GLONASS)
    source = 3;
end

% call prn2x function for GLONASS constellation
freqChan = navsu.svprn.prn2x(14, prn, epoch, 2*ones(size(prn)), source);

end
