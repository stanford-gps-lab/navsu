function [H,delta_z,residsPost,K,measMat,pred_meas,measMatRemoved,R] =  ...
    measUpdateExclude(H,cov_propagated,R,measMat,pred_meas,PARAMS)

% loop through and remove bad measurements, starting with the bad ones,
% until they are all clean

measMatRemoved = zeros(0,size(measMat,2));

largeResids  = true;
mediumResids = true;

idx = 1;
while largeResids || mediumResids
    % Keep iterating until there are no bad measurements
    
    % 7. Calculate Kalman gain using (3.21)
    K = cov_propagated * H' /(H *cov_propagated * H' + R);
    
    % 8. Measurement innovations
    fullMeas = measMat(:,5);
    delta_z  = fullMeas-pred_meas;
    
    residsPost = delta_z-H*K*delta_z;
    
    % Set the thresholds- remove large errors first
    %    if largeResids
    excludeThreshLarge = PARAMS.measUse.excludeThreshLarge(measMat(:,6))';
    %    else
    excludeThreshMedium = PARAMS.measUse.excludeThresh(measMat(:,6))';
    %    end
    
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
        measMatRemoved = [measMatRemoved; measMat(indsLargeResids,:)];
        
        % remove from H,R,measMat,pred_meas
        H(indsLargeResids,:) = [];
        R(indsLargeResids,:) = [];
        R(:,indsLargeResids) = [];
        measMat(indsLargeResids,:) = [];
        pred_meas(indsLargeResids,:) = [];
        
        
    elseif ~isempty(indsMediumResids) && largeResids == false
        % Only if we have found that there are no large residuals can we
        % remove using the small threshold.
         
        % save what we're removing
        measMatRemoved = [measMatRemoved; measMat(indsMediumResids,:)];
        
        % remove from H,R,measMat,pred_meas
        H(indsMediumResids,:) = [];
        R(indsMediumResids,:) = [];
        R(:,indsMediumResids) = [];
        measMat(indsMediumResids,:) = [];
        pred_meas(indsMediumResids,:) = [];
    end
    
    idx = idx+1;
    
    if idx > 20
        error('Measurement update is taking too long- please come check this out');
    end

end



end
