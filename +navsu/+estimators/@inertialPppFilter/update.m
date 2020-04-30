function update(obj,epoch,obs,corrData)


% there are essentially two kinds of updates that can happen here- IMU
% update and non-IMU update
% IMU updates only run the mechanization
% other measurement updates can run the entire thing

% If there is an IMU measurement in here, then
measTypes = cellfun(@(x) getfield(x,'type'),obs);

% run the inertial mechanization to propagate to the current time step
obj.mechanization(epoch,obs);

if length(measTypes) == 1 && ismember(measTypes,navsu.internal.MeasEnum.IMU)
    return;
end

% Manage the states in the filter :)
measRemovedSlip = navsu.ppp.manageStatesMulti(obj,epoch,obs);

% Time update (this is update of non-position, velocity, and attitude
% states as well as covariance propagation)
obj.timeUpdate(epoch)

% Measurement update
obj.measUpdate(epoch,obs,corrData,measRemovedSlip);

end

