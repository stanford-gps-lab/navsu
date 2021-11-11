function [measId,extraInputs] = update(obj, epoch, obs, corrData, varargin)

%%
% Manage the states in the filter :)
measRemovedSlip = obj.manageStatesMulti(epoch, obs);

% Time update
obj.timeUpdate(epoch)

% Measurement update
[measId, extraInputs] = obj.measUpdate(epoch, obs, corrData, measRemovedSlip, ...
                                       varargin{:});

% Make sure that the filter knows that it is running.
if (epoch - obj.epochLastGnssUpdate) < 10
    % if it's been too long, assume we're not initialized anymore
    obj.initialized = 2;
else
    obj.initialized = 0;
end

end
