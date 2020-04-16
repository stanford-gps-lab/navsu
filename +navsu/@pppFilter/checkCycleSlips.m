function measRemovedFull = checkCycleSlips(obj,epoch,gnssMeas,PARAMS)



% Compute current geometry free combination


measRemovedFull = [];

switch PARAMS.measUse.slipDetector
    case 'GFREE'
        % Check what carrier phase measurements are currently available
        indsCpDfAvail = find(gnssMeas.range.ind == 2 & gnssMeas.range.obs > 0 & gnssMeas.range.sig >= 100);
        
        % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsCpDfAvail) gnssMeas.range.constInds(indsCpDfAvail) ones(size(indsCpDfAvail)) gnssMeas.range.sig(indsCpDfAvail)];
        
        
        [~,indGFrees] = ismember(measInfoAvail,obj.cycleSlipInfo.measInfoGFree,'rows');
        
        % loop through and update each of the geometry free combinations
        for idx = 1:length(indsCpDfAvail)
            % Find each of the other measurements
            
            measSig1 = floor(measInfoAvail(idx,4)/100);
            measSig2 = measInfoAvail(idx,4)-measSig1*100;
            
            ind1 = find(gnssMeas.range.PRN == measInfoAvail(idx,1) & gnssMeas.range.constInds == measInfoAvail(idx,2) & ...
                gnssMeas.range.ind == 2 & gnssMeas.range.sig == measSig1);
            
            ind2 = find(gnssMeas.range.PRN == measInfoAvail(idx,1) & gnssMeas.range.constInds == measInfoAvail(idx,2) & ...
                gnssMeas.range.ind == 2 & gnssMeas.range.sig == measSig2);
            
            gFreei = gnssMeas.range.obs(ind1)-gnssMeas.range.obs(ind2);
            
            % Look for this in the geometry free combinations we already have
            indGFree = indGFrees(idx);
            
            if indGFree ~= 0
                % Check if we have a cycle slip
                gFreeLast = obj.cycleSlipInfo.gFree(indGFree);
                epochLast = obj.cycleSlipInfo.epochLastGFree(indGFree);
                
                cycleSlip = abs(gFreei-gFreeLast) > PARAMS.measUse.gFreeSlipThresh | (epoch-epochLast) > 30;
                
                if cycleSlip
                    % Just going to remove it from the state and covariance
                    % remove each of the
                    measInfoRemove = [measInfoAvail(idx,:); [measInfoAvail(idx,1:3) measSig1]; ...
                        measInfoAvail(idx,1:3) measSig2];
                    obj.removeFlexState(measInfoRemove)
                    
                    measRemovedFull = [measRemovedFull; measInfoRemove];
                end
                
                % Update the list
                obj.cycleSlipInfo.gFree(indGFree) = gFreei;
                obj.cycleSlipInfo.epochLastGFree(indGFree) = epoch;
                
            else
                % Add this one to the list
                obj.cycleSlipInfo.gFree          = [obj.cycleSlipInfo.gFree; gFreei];
                obj.cycleSlipInfo.epochLastGFree = [obj.cycleSlipInfo.epochLastGFree; epoch];
                obj.cycleSlipInfo.measInfoGFree  = [obj.cycleSlipInfo.measInfoGFree; measInfoAvail(idx,:)];
            end
        end
        
    case 'RX_OUTPUT'
        % Check if there are any slips
        %         dt = epoch-obj.cycleSlipInfo.epochLastGFree;
        % Check what carrier phase measurements are currently available
        indsCpAvail = find(gnssMeas.range.ind == 2 & gnssMeas.range.obs > 0);
        
        % PRN | CONST | STATE TYPE (1 = CP) | SIGNAL NUMBER
        measInfoAvail = [gnssMeas.range.PRN(indsCpAvail) gnssMeas.range.constInds(indsCpAvail) ones(size(indsCpAvail)) gnssMeas.range.sig(indsCpAvail)];
        
        if ~isempty(obj.cycleSlipInfo.gFree)
            tLockLast = obj.cycleSlipInfo.gFree(indsCpAvail);
            tLockCurr = gnssMeas.range.lockTime(indsCpAvail);
            
            indsSlip = find(tLockCurr < tLockLast | isnan(tLockLast));
            
            measRemovedFull = measInfoAvail(indsSlip,:);
        end
        
        for idx = 1:size(measRemovedFull,1)
            obj.removeFlexState(measRemovedFull(idx,:))
        end
        
%         else
        obj.cycleSlipInfo.epochLastGFree = epoch;
        obj.cycleSlipInfo.gFree = gnssMeas.range.lockTime;
end











end