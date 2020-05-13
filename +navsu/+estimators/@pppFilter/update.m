function measId = update(obj,epoch,obs,corrData)

% Manage the states in the filter :)
measRemovedSlip = navsu.ppp.manageStatesMulti(obj,epoch,obs);

% Time update
obj.timeUpdate(epoch)

% Measurement update
obj.measUpdate(epoch,obs,corrData,measRemovedSlip);

% Make sure that the filter knows that it is running.
obj.initialized = 2;

end






















