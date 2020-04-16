function removeFlexState(obj,measInfoRemove)

% NaNs are wildcards
% this
stateNotNeeded = find(ismember(obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),measInfoRemove,'rows'));


% 

if ~isempty(stateNotNeeded)
    % Remove from the state estimate, covariance, FLEX_STATES, and
    % FLEX_STATES_INFO
    obj.state(obj.INDS_STATE.FLEX_STATES(stateNotNeeded)) = [];
    obj.cov(obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
    obj.cov(:,obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
    obj.INDS_STATE.FLEX_STATES_INFO(stateNotNeeded,:) = [];
    obj.INDS_STATE.FLEX_STATES = obj.INDS_STATE.FLEX_STATE_MIN-1+...
        [1:size(obj.INDS_STATE.FLEX_STATES_INFO,1)]';
end





end