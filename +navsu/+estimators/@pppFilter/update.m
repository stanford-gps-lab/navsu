function measId = update(obj,epoch,obs,corrData,varargin)

p = inputParser;

p.addParameter('measExclude',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
measExclude = res.measExclude;
%%

% Manage the states in the filter :)
measRemovedSlip = navsu.ppp.manageStatesMulti(obj,epoch,obs);

% Time update
obj.timeUpdate(epoch)

% Measurement update
measId = obj.measUpdate(epoch,obs,corrData,measRemovedSlip,'measExclude',measExclude);

% Make sure that the filter knows that it is running.
obj.initialized = 2;

end






















