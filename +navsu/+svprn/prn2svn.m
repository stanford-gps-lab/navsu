function svn = prn2svn(prn, epoch, const, source)
% prn2svn
% DESCRIPTION:
%   Map from PRN to SVN at a given time.
%
% INPUT:
%   prn     - PRN
%   epoch   - time is in fractional year (e.g 2012.233 very crudely) or 
%             Julian day (e.g. 2444239.5)
%   
% OUTPUT:
%   svn     - satellite vehicle number.
%
% See also: navsu.svprn.svn2prn, navsu.time.cal2jd

if nargin < 2
    error('You must specify an svn and a time');
end

% Constellations
% GPS = 1, GLO = 2, GAL = 3, BDS = 4
if nargin < 3
    % Default to GPS
    const =  1;
end



% Check for string constellation type input and convert to number
if ischar(const)
    constNames = {'GPS','GLO','GAL','BDS','SBAS','QZSS'};
    const = find(~cellfun(@isempty,(strfind(constNames,const))));
    if isempty(const)
        warning('Invalid constellation name- defaulting to GPS')
        const = 1;
    end
end

if nargin < 4
    % Source of database (only applies to GLONASS)
    if const == 2
        source = 3;
    else
        source = 1;
    end
end

if any(epoch < 3000)
    sdx = epoch < 3000;
    year = floor(epoch(sdx));
    dayn = ceil(navsu.time.YearDays(year)*(epoch(sdx) - year));
    epoch(sdx) = navsu.time.doy2jd(year, dayn);
end
if length(epoch) == 1
    epoch = epoch*ones(size(prn));
end

% Pull svndata table
svndata = navsu.svprn.constSvnData(const,source);

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
        svn(idx) = svndata(sdx,1);
    else
        svn(idx) = NaN;
    end
end

end