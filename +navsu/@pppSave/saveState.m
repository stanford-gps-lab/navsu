function saveState(obj,user,PARAMS,varargin)

%% Parse inputs

p = inputParser;

p.addParameter('tdx',[]);
p.addParameter('epoch',[]);
p.addParameter('pl',[]);
p.addParameter('plLoc',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
tdx        = res.tdx;  % Index of full state saving matrices
epoch      = res.epoch; % Epoch to save
pl         = res.pl;
plLoc      = res.plLoc;

if isempty(tdx)
    error('do something about this please')
end

%% Save off the data

llhi = navsu.geo.xyz2llh(user.pos');

att = navsu.geo.dcm2euler123(user.R_b_e);

% compute ENU euler angles
[~,Ri] = navsu.geo.xyz2enu(user.vel',llhi(1)*pi/180,llhi(2)*pi/180);
%     attEnu = initializeAttitude(velEnu,PARAMS);
attEnu = navsu.geo.dcm2euler123(Ri*user.R_b_e);

covEnu = Ri*user.cov(obj.INDS_STATE.POS,obj.INDS_STATE.POS)*Ri';

nEpochsSmall = 1000;
if isempty(obj.stateSaveSmall)
    % Start a new small saving container- this will be continually filled
    % and dumped so as to not have to always deal with the large output
    % matrix
   obj.tdxSmall       = tdx:(tdx+nEpochsSmall-1); 
   obj.stateSaveSmall = nan(19,nEpochsSmall);
   obj.covSaveSmall   = nan(size(obj.covSaveFull,1),nEpochsSmall);
   obj.covEnuSmall    = nan(3,nEpochsSmall);
   obj.plSmall        = nan(3,nEpochsSmall);
   obj.plLocSmall     = nan(3,nEpochsSmall);
   
   obj.epochsSmall    = nan(1,nEpochsSmall);
   obj.tdx2           = 1;
end

% put the state into this new matrix
savei = [user.pos+user.R_b_e*PARAMS.IMU_ARM;  user.vel; att; attEnu; user.imuBiasStates; ...
    user.clockBias(1)/PARAMS.c;];

obj.epochsSmall(obj.tdx2) = obj.epochSave(tdx);
obj.stateSaveSmall(:,obj.tdx2) = savei;
obj.covSaveSmall(:,obj.tdx2) = diag(user.cov(1:(user.INDS_STATE.FLEX_STATE_MIN-1),1:(user.INDS_STATE.FLEX_STATE_MIN-1)));
obj.tdxSmall(obj.tdx2) = tdx;
obj.covEnuSmall(:,obj.tdx2) = diag(covEnu);
if ~isempty(pl)
    obj.plSmall(:,obj.tdx2) = pl;
end
if ~isempty(plLoc)
    obj.plLocSmall(:,obj.tdx2) = plLoc;
end

if obj.tdx2 == nEpochsSmall
    % we're full- time to dump
    obj.stateSaveFull(:,obj.tdxSmall) = obj.stateSaveSmall;
    obj.covSaveFull(:,obj.tdxSmall)   = obj.covSaveSmall;
    obj.epochsFull(obj.tdxSmall)      = obj.epochsSmall;
    obj.covEnuFull(:,obj.tdxSmall)    = obj.covEnuSmall;
    obj.plFull(:,obj.tdxSmall)        = obj.plSmall;
    obj.plLocFull(:,obj.tdxSmall)     = obj.plLocSmall;
    
    obj.tdxSmall = [];
    obj.stateSaveSmall = [];
    obj.epochsSmall = [];
    obj.tdx2 = [];
    obj.covEnuSmall = [];
    obj.plSmall = [];
    obj.plLocSmall = [];
else
    obj.tdx2 = obj.tdx2+1;
end
    
    

% Save covariance
% lol = diag(user.cov(1:(user.INDS_STATE.FLEX_STATE_MIN-1),1:(user.INDS_STATE.FLEX_STATE_MIN-1)));
% obj.covSave(tdx,:) = lol;


end