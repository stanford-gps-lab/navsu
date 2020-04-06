function initDcb(obj,years,doys,varargin)

p = inputParser;
p.addParameter('source',  []);  % DLR
p.addParameter('FLAG_NO_LOAD', false);

settings = obj.settings;
% parse the results
parse(p, varargin{:});

res = p.Results;

if isempty(res.source)
    source = settings.dcbSource;
else
    source = res.source;
end
FLAG_NO_LOAD = res.FLAG_NO_LOAD;

% if there is no precise ephemeris already there, just make a new one
if isempty(obj.dcb)
    dcb = utility.readfiles.loadDcb(years,doys,settings,FLAG_NO_LOAD,source);
    obj.dcb = dcb;
else
    % This will be added to the rest of it
    
    % ADD STUFF HERE PLEASE
end




end