function update(obj,epoch,obs,corrData)


% there are essentially two kinds of updates that can happen here- IMU
% update and non-IMU update
% IMU updates only run the mechanization
% other measurement updates can run the entire thing

% Manage the states in the filter :)
% measRemovedSlip = navsu.ppp.manageStatesMulti(obj,epoch,obs);

% Time update
% obj.timeUpdate(epoch)

% Measurement update
% obj.measUpdate(epoch,obs,corrData,measRemovedSlip);

end

