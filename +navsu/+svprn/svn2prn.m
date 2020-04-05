function prn = svn2prn(svn, epoch,const)
%SVN = SVN2PRN(SVN, EPOCH)
%converts svn number to prn numbers given the time
% time is in fractional year (e.g 2012.233 very crudely) or
% Julian day (e.g. 2444239.5)
% CONST is index of GNSS constellation:
% GPS = 1, GLO = 2, GAL = 3, BDS = 4
%   If not included, defaults to GPS
%SEE ALSO CAL2JD, DOY2JD PRN2SVN


if nargin < 2
    error('You must specify a prn and a time');
end

% Constellations
% GPS = 1, GLO = 2, GAL = 3, BDS = 4
if nargin < 3
   const =  1;
end
% Check for string constellation type input and convert to number
if ischar(const)
    constNames = {'GPS','GLO','GAL','BDS'};
    const = find(~cellfun(@isempty,(strfind(constNames,const))));
    if isempty(const)
        warning('Invalid constellation name- defaulting to GPS')
        const = 1;
    end
end


if any(epoch < 3000)
    pdx = epoch < 3000;
    year = floor(epoch(pdx));
    dayn = ceil(navsu.time.yearDays(year)*(epoch(pdx) - year));
    epoch(pdx) = navsu.time.doy2jd(year, dayn);
end
if length(epoch) == 1
    epoch = epoch*ones(size(svn));
end

% Pull svndata table 
svndata = navsu.svprn.constSvnData(const);

%find start date in jd       
svndata(:,end+1) = navsu.time.cal2jd(svndata(:,3), svndata(:,4), ...
                    svndata(:,5) + (svndata(:,6) + svndata(:,7)/60)/24);
%find end date in jd       
svndata(:,end+1) = navsu.time.cal2jd(svndata(:,8), svndata(:,9), ...
                    svndata(:,10) + (svndata(:,11) + svndata(:,12)/60)/24);       
% fix infinities
svndata(svndata(:,8) == Inf,end) = Inf;
                
prn = NaN(size(svn));       

for sdx = 1:size(svndata,1)
    pdx = svndata(sdx,1) == svn & ...
           epoch >= svndata(sdx,end-1) & epoch <= svndata(sdx,end); 
    if any(pdx)
        prn(pdx) = svndata(sdx,2);
    end
end

end