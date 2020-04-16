function saveMeasRemoved(obj,epoch,measLow,measResids,measSlip)

if nargin < 5
    measSlip = zeros(0,4);
end

nMeasPreallocate = 1e3;

if isempty(obj.indRemoved)
   % Initialize
   obj.indRemoved    = 1;
   obj.measRemoved   = nan(nMeasPreallocate,5);
   obj.epochsRemoved = nan(nMeasPreallocate,1);
end
    
% Actually only keep the prn, const, and reason for elevation removals
measLow = unique(measLow(:,[1 2]),'rows');

% PRN | CONST | SIG  | MEAS TYPE (1=CODE,2=CARR,3=DOP) | REMOVAL REASON (1=ELEVATION,2=RESIDUALS, 3 = SLIP)
measRemove = [measLow(:,1:2) 0*ones(size(measLow,1),2) 1*ones(size(measLow,1),1);
             measResids(:,[1:3 6]) 2*ones(size(measResids,1),1);
             measSlip(:,[1 2 4]) 2*ones(size(measSlip,1),1) 3*ones(size(measSlip,1),1)];

% Check if we have space to add these to the matrix
nMeasRemove = size(measRemove,1);

if nMeasRemove+obj.indRemoved >= size(obj.measRemoved,1)
    % Need to add more space
   obj.measRemoved = [obj.measRemoved; nan(nMeasPreallocate,5)];
   obj.epochsRemoved = [obj.epochsRemoved; nan(nMeasPreallocate,1)];
end

% Add the newly removed measurements
indsAdd = obj.indRemoved:(obj.indRemoved+nMeasRemove-1);
obj.measRemoved(indsAdd,:) = measRemove;
obj.epochsRemoved(indsAdd) = epoch;

obj.indRemoved = obj.indRemoved+nMeasRemove;



end