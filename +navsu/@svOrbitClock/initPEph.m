function initPEph(obj,years,doys,varargin)

p = inputParser;

p.addParameter('FLAG_NO_LOAD',false);
p.addParameter('atxData',[]);
p.addParameter('FLAG_APC_OFFSET',false);
p.addParameter('TIME_STRIP',true);

settings = obj.settings;

% parse the results
parse(p, varargin{:});
res = p.Results;
FLAG_NO_LOAD = res.FLAG_NO_LOAD;       % Flag to not actually load the data- should probably not be used here
atxData      = res.atxData;            % IGS ANTEX data
FLAG_APC_OFFSET = res.FLAG_APC_OFFSET; % whether or not to add the antenna phase center offset
TIME_STRIP      = res.TIME_STRIP;      % strip output data to only current day

% if there is no precise ephemeris already there, just make a new one
if isempty(obj.PEph)
    Peph = utility.readfiles.loadPEph(years,doys,settings,FLAG_NO_LOAD,atxData,FLAG_APC_OFFSET,TIME_STRIP);
    
    prnConstInds = unique([Peph.PRN Peph.constellation],'rows');
    prns = prnConstInds(:,1);
    constInds = prnConstInds(:,2);
    
    epochsPeph = unique(Peph.epochs);
    Peph2.clock_bias    = zeros(length(epochsPeph),length(prns));
    Peph2.clock_drift   = zeros(length(epochsPeph),length(prns));
    Peph2.position      = zeros(length(epochsPeph),3,length(prns));
    Peph2.velocity      = zeros(length(epochsPeph),3,length(prns));
    Peph2.event         = zeros(length(epochsPeph),length(prns));
    Peph2.epochs        = epochsPeph;
    Peph2.PRN           = prns;
    Peph2.constellation = constInds;
    for pdx = 1:length(prns)
        indsi = find(Peph.PRN == prns(pdx) & Peph.constellation == constInds(pdx));
        [~,inds2] = ismember(epochsPeph,Peph.epochs(indsi));
        
        Peph2.clock_bias(inds2,pdx) = Peph.clock_bias(indsi);
        Peph2.clock_drift(inds2,pdx) = Peph.clock_drift(indsi);
        Peph2.position(inds2,:,pdx) = Peph.position(indsi,:);
        Peph2.velocity(inds2,:,pdx) = Peph.velocity(indsi,:);
        Peph2.Event(inds2,pdx)      = Peph.Event(indsi);
    end
    
    obj.PEph = Peph2;

else
    % This will be added to the rest of it
    
    % ADD STUFF HERE PLEASE
end




end