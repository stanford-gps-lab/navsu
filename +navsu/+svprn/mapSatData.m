function svn = mapSatData(queryCol, refCol, satId, epoch, const, source)
% getSatData(queryCol, refCol, satId, epoch, const, source)
% DESCRIPTION:
%   Map from one satellite specific value toanother at a given time.
%
% INPUT:
%   queryCol- index of column of desired data from the svn table (see
%             navsu.svprn.constSvnData)
%   refCol  - index of column of reference data from the svn table (see
%             navsu.svprn.constSvnData)
%   satId   - value of the reference data (e.g. PRN, SVN, ...)
%   epoch   - time is in fractional year (e.g 2012.233 very crudely) or 
%             Julian day (e.g. 2444239.5)
%   const   - N x 1 vector of constellation indices
%   source  - OPTIONAL data source (only applies to GLONASS)
%   
% OUTPUT:
%   svn     - satellite vehicle number.
%
% See also: navsu.svprn.prn2svn, navsu.svprn.prn2FreqChanGlonass

if nargin < 4
    error('You must specify a prn and a time');
end

if isrow(satId)
    % ensure we have a column vector of PRNs
    satId = satId';
end

multiEpoch = numel(epoch) > 1;

if multiEpoch && isrow(epoch)
    % need scalar or column vector
    epoch = epoch';
end

% Constellations
% GPS = 1, GLO = 2, GAL = 3, BDS = 4
if nargin < 5
    % Default to GPS
    const = ones(size(satId));
elseif ischar(const)
    % Check for string constellation type input and convert to number
    const = navsu.svprn.convertConstIndName(const);
    if isempty(const)
        warning('Invalid constellation name - defaulting to GPS')
        const = 1;
    end
elseif isrow(const)
    const = const';
end

if numel(const) == 1
    const = const * ones(size(satId));
end

if nargin < 6
    % Source of database (only applies to GLONASS)
    source = ones(size(const));
    source(const == 2) = 3;
elseif numel(source) == 1
    source = source * ones(size(const));
elseif isrow(source)
    source = source';
end

% convert fractional year to julian date if necessary
sdx = epoch < 3000;
if any(sdx)
    year = floor(epoch(sdx));
    dayn = ceil(navsu.time.YearDays(year)*(epoch(sdx) - year));
    epoch(sdx) = navsu.time.doy2jd(year, dayn);
end


svn = NaN(size(satId));

% get constellation - source combinations to consider
unCS = unique([const, source], 'rows');

for csdx = 1:size(unCS, 1)
    % retrieve constellation and source
    consti = unCS(csdx, 1);
    sourcei = unCS(csdx, 2);
    
    % Pull corresponding svndata table
    [svndata, blockText] = navsu.svprn.constSvnData(consti, sourcei);
    
    % find start date in jd
    epStart = navsu.time.cal2jd_vect(svndata(:,3), svndata(:,4), ...
        svndata(:,5)) + (svndata(:,6) + svndata(:,7)/60)/24;
    % find end date in jd
    epEnd = navsu.time.cal2jd_vect(svndata(:,8), svndata(:,9), ...
        svndata(:,10)) + (svndata(:,11) + svndata(:,12)/60)/24;
    % fix infinities
    epEnd(svndata(:,8) == Inf) = Inf;
    
    % save data for satellites that are part of this case
    
    idx = find(const == consti & source == sourcei);
    sati = satId(idx)';
    if multiEpoch
        epochi = epoch(idx)';
    else
        epochi = epoch;
    end
    sdx = sati == svndata(:, refCol) ...
        & epochi >= epStart ...
        & epochi < epEnd;
    % which PRNs do I have data for?
    haveData = any(sdx, 1);
    % get indices in SVN table, store values from desired column
    sdxInd = find(sdx) - size(svndata, 1)*(find(haveData)'-1);
    svn(idx(haveData)) = svndata(sdxInd, queryCol);
    
end

end