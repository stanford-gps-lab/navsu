function initPClock(obj,years,doys,varargin)

p = inputParser;
% 
% p.addParameter('FLAG_NO_LOAD',false);
% p.addParameter('atxData',[]);
% p.addParameter('FLAG_APC_OFFSET',true);

settings = obj.settings;
% parse the results
parse(p, varargin{:});
res = p.Results;
% FLAG_NO_LOAD = res.FLAG_NO_LOAD;       % Flag to not actually load the data- should probably not be used here
% atxData      = res.atxData;            % IGS ANTEX data
% FLAG_APC_OFFSET = res.FLAG_APC_OFFSET; % whether or not to add the antenna phase center offset

% if there is no precise ephemeris already there, just make a new one
if isempty(obj.PClock)
    PClock = utility.readfiles.loadCFst(years,doys,settings);
    obj.PClock = PClock;
else
    % This will be added to the rest of it
    
    % ADD STUFF HERE PLEASE
end




end