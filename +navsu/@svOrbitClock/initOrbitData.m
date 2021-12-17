function initOrbitData(obj,years,doys,varargin)
% Download (if necessary) and load GNSS orbit data
% DESCRIPTION:
%   Initialize the precise orbit structure from navsu.svOrbitClock.  This
%   will load from local files or download from IGS FTP sites if necessary.
%   
% INPUT:
%   years   - vector of years associated with the days of years when data
%             should be downloaded
%   doys    - vector of days of year for which data should be downloaded
%  
% OPTIONAL INPUTS:
%   atxData - antenna phase center structure parsed from IGS .atx file-
%             this would be the output of navsu.svOrbitClock.initAtxData.
%             Used if you want to offset the satellite positions from
%             center of mass to antenna phase center right now.  You
%             probably can do this later rather than now!
%   FLAG_APC_OFFSET - boolean indicating whether or not to offset the
%             satellite centers of mass to their antenna phase centers.
%             Requires atxData.  Default is FALSE.
%   DOWNLOAD - boolean indcating whether or not to try to download the
%             required products if they aren't found locally. Default is
%             TRUE. 
%             
% OUTPUT:
%   The object will have a brand new orbit structure!  
%
% See also: navsu.svOrbitClock.propagate

%%
p = inputParser;
p.addParameter('atxData',[]);
p.addParameter('FLAG_APC_OFFSET',false);
p.addParameter('DOWNLOAD',true);   

% parse the results
parse(p, varargin{:});
res = p.Results;
atxData         = res.atxData;            % IGS ANTEX data
FLAG_APC_OFFSET = res.FLAG_APC_OFFSET; % whether or not to add the antenna phase center offset
DOWNLOAD        = res.DOWNLOAD;        % indicator to check for downloads and download

settings = obj.settings;

FLAG_NO_LOAD = false;       % Flag to not actually load the data- should probably not be used here

% if there is no precise ephemeris already there, just make a new one
if isempty(obj.PEph)
    
    if DOWNLOAD
        [~,filenames,filenameFull] = navsu.readfiles.loadPEph( ...
            years, doys, settings, true, atxData, FLAG_APC_OFFSET);
        
        % build a list of IGS AC codes to download
        fileAvailable = cellfun(@exist,filenameFull);
        igsCodes = cellfun(@(x) x(1:3),filenames,'UniformOutput',false);
        
        codesDownload = unique(igsCodes(~fileAvailable));
        
        for idx =1 :length(codesDownload)
           navsu.ftp.download(1,years,doys,settings,codesDownload{idx});
        end
    end
    
    Peph = navsu.readfiles.loadPEph(years,doys,settings,FLAG_NO_LOAD,atxData,FLAG_APC_OFFSET);
    
    % exclude values where the PRN is NaN
    nPRNvals = length(Peph.PRN);
    finPRN = isfinite(Peph.PRN);

    fn = fieldnames(Peph);
    for fni = 1:length(fn)

        if size(Peph.(fn{fni}), 1) == nPRNvals
            Peph.(fn{fni}) = Peph.(fn{fni})(finPRN, :);
        end

    end
    
    prnConstInds = unique([Peph.PRN Peph.constellation], 'rows');
    prns = prnConstInds(:,1);
    constInds = prnConstInds(:,2);
    
    epochsPeph = unique(Peph.epochs);
    Peph2.clock_bias  = NaN(length(epochsPeph), length(prns));
    Peph2.clock_drift = NaN(length(epochsPeph), length(prns));
    Peph2.position    = NaN(length(epochsPeph), 3, length(prns));
    Peph2.velocity    = NaN(length(epochsPeph), 3, length(prns));
    Peph2.event       = NaN(length(epochsPeph), length(prns));
    Peph2.epochs      = epochsPeph;
    Peph2.PRN = prns;
    Peph2.constellation = constInds;
    for pdx = 1:length(prns)
        indsi = Peph.PRN == prns(pdx) & Peph.constellation == constInds(pdx);
%         [~,inds2] = ismember(epochsPeph,Peph.epochs(indsi));
        [~,inds2] = ismember(Peph.epochs(indsi), epochsPeph);
        
        Peph2.clock_bias(inds2, pdx) = Peph.clock_bias(indsi);
        Peph2.clock_drift(inds2, pdx) = Peph.clock_drift(indsi);
        Peph2.position(inds2, :, pdx) = Peph.position(indsi,:);
        Peph2.velocity(inds2, :, pdx) = Peph.velocity(indsi,:);
        Peph2.Event(inds2, pdx)      = Peph.Event(indsi);
    end
    
    obj.PEph = Peph2;

else
    % This will be added to the rest of it
    
    % ADD STUFF HERE PLEASE
end




end