function svn = prn2FreqChanGlonass(prn, epoch, source)
%PRN = PRN2SVN(PRN, EPOCH,CONST)
%converts prn number to svn numbers given the time
% time is in fractional year (e.g 2012.233 very crudely) or
% Julian day (e.g. 2444239.5)
% CONST is index of GNSS constellation:
% GPS = 1, GLO = 2, GAL = 3, BDS = 4
%   If not included, defaults to GPS
%SEE ALSO CAL2JD, DOY2JD SVN2PRN

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
                
svn = NaN(size(prn));       

% for pdx = 1:size(svndata,1)
%     sdx = svndata(pdx,2) == prn & ...
%            epoch >= svndata(pdx,end-1) & epoch <= svndata(pdx,end); 
%     if any(sdx)
%         svn(sdx) = svndata(pdx,1);
%         svn(sdx) = length(pdx);
%     end
% end

for idx = 1:length(prn)

    prni = prn(idx);
    epochi = epoch(idx);
    sdx = find(svndata(:,2) == prni & epochi >= svndata(:,end-1) & epochi < svndata(:,end));
    if ~isempty(sdx)
        svn(idx) = svndata(sdx,14);
    else
        svn(idx) = NaN;
    end
end

end