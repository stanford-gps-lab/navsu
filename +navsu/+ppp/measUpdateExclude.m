function [H,delta_z,residsPost,K,predMeas,measMatRemoved,R,measIdRemoved,measId] =  ...
    measUpdateExclude(H,cov_propagated,R,predMeas,PARAMS,meas,measId)

% loop through and remove bad measurements, starting with the bad ones,
% until they are all clean

measMatRemoved = zeros(0,6);
measIdRemoved = [];

largeResids  = true;
mediumResids = true;

% Build exclusion factor vector
threshLarge = zeros(size(measId, 1), 1);
threshMedium = zeros(size(measId, 1), 1);
measTypeList = cat(1, measId.TypeID);
types = unique(measTypeList);

for idx = 1:length(types)
    if types(idx) == navsu.internal.MeasEnum.GNSS
        indsGnss = find(types(idx) == measTypeList);
        
        measSubtypes = cat(1,measId(indsGnss).subtype);
        indsCode = indsGnss(measSubtypes == navsu.internal.MeasEnum.Code);
        indsCarrier = indsGnss(measSubtypes == navsu.internal.MeasEnum.Carrier);
        indsDoppler = indsGnss(measSubtypes == navsu.internal.MeasEnum.Doppler);

        threshLarge(indsCode) = PARAMS.measUse.excludeThreshLarge.GNSS.Code;
        threshLarge(indsCarrier) = PARAMS.measUse.excludeThreshLarge.GNSS.Carrier;
        threshLarge(indsDoppler) = PARAMS.measUse.excludeThreshLarge.GNSS.Doppler;
        
        threshMedium(indsCode) = PARAMS.measUse.excludeThresh.GNSS.Code;
        threshMedium(indsCarrier) = PARAMS.measUse.excludeThresh.GNSS.Carrier;
        threshMedium(indsDoppler) = PARAMS.measUse.excludeThresh.GNSS.Doppler;
    else
        threshLarge(types(idx) == measTypeList) = PARAMS.measUse.excludeThreshLarge.(char(types(idx)));
        threshMedium(types(idx) == measTypeList) = PARAMS.measUse.excludeThresh.(char(types(idx)));
    end
end



idx = 1;
while largeResids || mediumResids
    % Keep iterating until there are no bad measurements
    
    % 7. Calculate Kalman gain using (3.21)
    K = (cov_propagated * H') /(H *cov_propagated * H' + R);
    
    % 8. Measurement innovations
    delta_z  = meas - predMeas;
    
    residsPost = delta_z - H*K*delta_z;
    
    % Set the thresholds- remove large errors first
     
    indsLargeResids  = find(abs(residsPost) > threshLarge);
    indsMediumResids = find(abs(residsPost) > threshMedium);
    
    % Check where we are currently- are there large or small residuals?
    largeResids = ~isempty(indsLargeResids);
    mediumResids = ~isempty(indsMediumResids);
    
    % can only remove using EITHER large or small thresholds, not one then
    % the other
    if largeResids
        % Remove large residuals
        indsRemove = indsLargeResids;
        
    elseif mediumResids
        % Only if we have found that there are no large residuals can we
        % remove using the small threshold.
        indsRemove = indsMediumResids;
        
    else
        indsRemove = [];
    end
    % save what we're removing
    measIdRemoved = [measIdRemoved; measId(indsRemove)];
        
    % remove from H, R, measMat, pred_meas
    H(indsRemove, :) = [];
    R(indsRemove, :) = [];
    R(:, indsRemove) = [];
    predMeas(indsRemove, :) = [];

    meas(indsRemove) = [];
    measId(indsRemove) = [];
    threshLarge(indsRemove) = [];
    threshMedium(indsRemove) = [];    
    
    idx = idx+1;
    
    if idx > 20
        error('Measurement update is taking too long- please come check this out');
    end

end



end
