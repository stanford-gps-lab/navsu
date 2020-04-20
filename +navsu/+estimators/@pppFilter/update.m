function update(obj,epoch,obs,corrData)

% Manage the states in the filter :)
measRemovedSlip = navsu.ppp.manageStatesMulti(obj,epoch,obs);

% Time update
obj.timeUpdate(epoch)

% Measurement update
obj.measUpdate(epoch,obs,corrData,measRemovedSlip);


end






















