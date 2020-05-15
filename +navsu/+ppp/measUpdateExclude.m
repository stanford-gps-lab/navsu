function [H,delta_z,residsPost,K,predMeas,measMatRemoved,R,measIdRemoved,measId] =  ...
    measUpdateExclude(H,cov_propagated,R,predMeas,PARAMS,meas,measId)

% loop through and remove bad measurements, starting with the bad ones,
% until they are all clean

measMatRemoved = zeros(0,6);
measIdRemoved = [];

largeResids  = true;
mediumResids = true;

idx = 1;
while largeResids || mediumResids
    % Keep iterating until there are no bad measurements
    
    % 7. Calculate Kalman gain using (3.21)
    K = cov_propagated * H' /(H *cov_propagated * H' + R);
    
    % 8. Measurement innovations
    fullMeas = meas;
    delta_z  = fullMeas-predMeas;
    
    residsPost = delta_z-H*K*delta_z;
    
    % Set the thresholds- remove large errors first
    excludeThreshLarge = PARAMS.measUse.excludeThreshLarge*sqrt(diag(R));
    excludeThreshMedium = PARAMS.measUse.excludeThresh*sqrt(diag(R));
%     
    indsLargeResids  = find(abs(residsPost)>excludeThreshLarge);
    indsMediumResids = find(abs(residsPost)>excludeThreshMedium);
    
    % Check where we are currently- are there large or small residuals?
    if isempty(indsLargeResids)
        largeResids = false;
    end
    
    if isempty(indsMediumResids)
        % If we have no residuals exceeding our smaller threshold, just move
        % forward with what we have- these are clean
        mediumResids = false;
    end
    
    % can only remove using EITHER large or small thresholds, not one then
    % the other
    if ~isempty(indsLargeResids)
        % Remove large residuals
        % save what we're removing
        
        measIdRemoved = [measIdRemoved; measId(indsLargeResids)];
        
        % remove from H,R,measMat,pred_meas
        H(indsLargeResids,:) = [];
        R(indsLargeResids,:) = [];
        R(:,indsLargeResids) = [];
        predMeas(indsLargeResids,:) = [];
        
        meas(indsLargeResids) = [];
        measId(indsLargeResids) = [];
        
        
    elseif ~isempty(indsMediumResids) && largeResids == false
        % Only if we have found that there are no large residuals can we
        % remove using the small threshold.
         
        % save what we're removing
        measIdRemoved = [measIdRemoved; measId(indsMediumResids)];
        
        % remove from H,R,measMat,pred_meas
        H(indsMediumResids,:) = [];
        R(indsMediumResids,:) = [];
        R(:,indsMediumResids) = [];
        predMeas(indsMediumResids,:) = [];
        
        meas(indsMediumResids) = [];
        measId(indsMediumResids) = [];
    end
    
    idx = idx+1;
    
    if idx > 20
        error('Measurement update is taking too long- please come check this out');
    end

end



end
