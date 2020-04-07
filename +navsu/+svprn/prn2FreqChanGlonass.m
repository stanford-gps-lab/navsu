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

if any(epoch < 3000)
    sdx = epoch < 3000;
    year = floor(epoch(sdx));
    dayn = ceil(navsu.time.yearDays(year)*(epoch(sdx) - year));
    epoch(sdx) = navsu.time.doy2jd(year, dayn);
end
if length(epoch) == 1
    epoch = epoch*ones(size(prn));
end

% Pull svndata table 
svndata = navsu.svprn.constSvnData(2,source);

%find start date in jd       
svndata(:,end+1) = navsu.time.cal2jd_vect(svndata(:,3), svndata(:,4), ...
                    svndata(:,5)) + (svndata(:,6) + svndata(:,7)/60)/24;
%find end date in jd       
svndata(:,end+1) = navsu.time.cal2jd_vect(svndata(:,8), svndata(:,9), ...
                    svndata(:,10) + (svndata(:,11) + svndata(:,12)/60)/24);       
% fix infinities
svndata(svndata(:,8) == Inf,end) = Inf;
                
freqChan = NaN(size(prn));       


for idx = 1:length(prn)
    prni = prn(idx);
    epochi = epoch(idx);
    sdx = find(svndata(:,2) == prni & epochi >= svndata(:,end-1) & epochi < svndata(:,end));
    if ~isempty(sdx)
        freqChan(idx) = svndata(sdx,14);
    else
        freqChan(idx) = NaN;
    end
end

end