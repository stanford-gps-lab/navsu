function closeSaveState(obj)

% Put the interim save matrix into the big one
if ~isempty(obj.tdx2)
    obj.stateSaveFull(:,obj.tdxSmall(1:(obj.tdx2-1))) = obj.stateSaveSmall(:,1:(obj.tdx2-1));
    obj.covSaveFull(:,obj.tdxSmall(1:(obj.tdx2-1))) = obj.covSaveSmall(:,1:(obj.tdx2-1));
    
    obj.epochsFull(obj.tdxSmall(1:(obj.tdx2-1))) = obj.epochsSmall(1:(obj.tdx2-1));
    obj.covEnuFull(:,obj.tdxSmall(1:(obj.tdx2-1))) = obj.covEnuSmall(:,1:(obj.tdx2-1));
    obj.plFull(:,obj.tdxSmall(1:(obj.tdx2-1))) = obj.plSmall(:,1:(obj.tdx2-1));
end

% Pull everything out
obj.posSave = obj.stateSaveFull(1:3,:)';
obj.velSave = obj.stateSaveFull(4:6,:)';
obj.attSave = obj.stateSaveFull(7:9,:)';
obj.attEnuSave = obj.stateSaveFull(10:12,:)';
obj.biasSave = obj.stateSaveFull(13:18,:)';
obj.clockSave = obj.stateSaveFull(19:end,:)';

obj.covSave  = obj.covSaveFull';
end