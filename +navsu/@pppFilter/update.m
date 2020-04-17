function [measMatRemoved,measMatRemovedLow] = ...
    update(obj,epoch,obs,corrData,outStruc)


R_b_e = obj.R_b_e;
pos   = obj.pos;
vel   = obj.vel;
imuBiasStates = obj.imuBiasStates;
cov  = obj.cov;
nState = size(cov,1);

c = navsu.constants.c;

PARAMS = obj.PARAMS;

if isempty(obj.lastAccMeas)
    accMeas = zeros(3,1);
else
    accMeas = obj.lastAccMeas;
end
if isempty(obj.lastGyroMeas)
    gyroMeas = zeros(3,1);
else
    gyroMeas = obj.lastGyroMeas;
end

llhi = navsu.geo.xyz2llh(pos');
latRad = llhi(1)*pi/180;

pos0 = obj.posPrevTc;

epoch0 = epoch;

% Time since last GNSS kalman filter update
dtKf = epoch-obj.epochLastGnssUpdate;
% Time since last inertial update
tSinceInertialUpdate = epoch-obj.epochLastInertialUpdate;

obj.epochLastGnssUpdate = epoch;

% Constants (sone of these could be changed to inputs at a later date)
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
R_0 = 6378137; %WGS84 Equatorial radius in meters
e = 0.0818191908425; %WGS84 eccentricity
if tSinceInertialUpdate < dtKf
    % Skew symmetric matrix of Earth rate
    Omega_ie = crossProdMatrix([0,0,omega_ie]);
    
    % 1. Determine transition matrix using (14.50) (first-order approx)
    Phi_matrix = eye(size(cov));
    
    Phi_matrix(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE) = ...
        Phi_matrix(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE) - Omega_ie * dtKf;
    
    Phi_matrix(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.W_BIAS) = ...
        R_b_e * dtKf;
    
    Phi_matrix(obj.INDS_STATE.VEL,obj.INDS_STATE.ATTITUDE) = ...
        -dtKf * crossProdMatrix(R_b_e * accMeas);
    
    Phi_matrix(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL) = ...
        Phi_matrix(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL) - 2 * Omega_ie * dtKf;
    
    geocentric_radius = R_0 / sqrt(1 - (e * sin(latRad))^2) *...
        sqrt(cos(latRad)^2 + (1 - e^2)^2 * sin(latRad)^2); % from (2.137)
    
    Phi_matrix(obj.INDS_STATE.VEL,obj.INDS_STATE.POS) = ...
        -dtKf * 2 * Gravity_ECEF(pos) /geocentric_radius*pos'/sqrt(pos'*pos);
    
    Phi_matrix(obj.INDS_STATE.VEL,obj.INDS_STATE.ACC_BIAS) = R_b_e * dtKf;
    
    Phi_matrix(obj.INDS_STATE.POS,obj.INDS_STATE.VEL) = eye(3) * dtKf;
    
    Phi_matrix(obj.INDS_STATE.CLOCK_BIAS,obj.INDS_STATE.CLOCK_DRIFT) = dtKf;
    
    % 2. Determine approximate system noise covariance matrix using (14.82)
    Q              = zeros(size(cov));
    
    if 1
        Q(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE) = ...
            eye(3) * PARAMS.Q.gyro_noise_PSD.^2 * dtKf;
        
        Q(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL)     =  ...
            eye(3) * PARAMS.Q.accel_noise_PSD.^2 * dtKf;
        
        Q(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.ACC_BIAS) = ...
            eye(3) * PARAMS.Q.ACC_BIAS.^2 * dtKf;
        
        Q(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.W_BIAS) = ...
            eye(3) * PARAMS.Q.W_BIAS.^2 * dtKf;
        
        Q(obj.INDS_STATE.CLOCK_BIAS,obj.INDS_STATE.CLOCK_BIAS) = ...
            eye(length(obj.INDS_STATE.CLOCK_BIAS)).*PARAMS.Q.RXB.^2 * dtKf/5000;
        
        Q(obj.INDS_STATE.CLOCK_DRIFT,obj.INDS_STATE.CLOCK_DRIFT) = ...
            eye(length(obj.INDS_STATE.CLOCK_BIAS)).*PARAMS.Q.RXB.^2 * dtKf;
    else
        % Don't approximate so much
        Srg = PARAMS.Q.gyro_noise_PSD.^2*0.01;
        Sra = PARAMS.Q.accel_noise_PSD.^2*0.01;
        
        Sbad = PARAMS.Q.ACC_BIAS.^2;
        Sbgd = PARAMS.Q.W_BIAS.^2;
        
        Ts = dtKf;
        
        F21 = -crossProdMatrix(R_b_e * accMeas);
        
        Q11 = (Srg*Ts+1/3*Sbgd*Ts^3)*eye(3);
        Q21 = (1/2*Srg*Ts^2+1/4*Sbgd*Ts^4)*F21;
        Q22 = (Sra*Ts+1/3*Sbad*Ts^3)*eye(3)+(1/3*Srg*Ts^3+1/5*Sbgd*Ts^5)*F21*F21'; %#ok<MHERM>
        Q31 = (1/3*Srg*Ts^3+1/5*Sbgd*Ts^5)*F21;
        Q32 = (1/2*Sra*Ts^2+1/4*Sbad*Ts^4)*eye(3)+(1/4*Srg*Ts^4+1/6*Sbgd*Ts^6)*F21*F21';
        Q33 = (1/3*Sra*Ts^3+1/5*Sbad*Ts^5)*eye(3)+(1/5*Srg*Ts^5+1/7*Sbgd*Ts^7)*F21*F21';
        Q34 = 1/3*Sbad*Ts^3*R_b_e;
        Q35 = 1/4*Sbgd*Ts^4*F21*R_b_e;
        
        Q15 = 1/2*Sbgd*Ts^2*R_b_e;
        Q24 = 1/2*Sbad*Ts^2*R_b_e;
        Q25 = 1/3*Sbgd*Ts^3*F21*R_b_e;
        Q44 = Sbad*Ts*eye(3);
        Q55 = Sbgd*Ts*eye(3);
        
        % 1 = attitude, 2 = velocity, 3 = position, 4 = acc bias, 5 = w bias
        Q(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE) = Q11;
        Q(obj.INDS_STATE.VEL,obj.INDS_STATE.ATTITUDE) = Q21;
        Q(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL) = Q22;
        Q(obj.INDS_STATE.POS,obj.INDS_STATE.ATTITUDE) = Q31;
        Q(obj.INDS_STATE.POS,obj.INDS_STATE.VEL) = Q32;
        Q(obj.INDS_STATE.POS,obj.INDS_STATE.POS) = Q33;
        Q(obj.INDS_STATE.POS,obj.INDS_STATE.ACC_BIAS) = Q34;
        Q(obj.INDS_STATE.POS,obj.INDS_STATE.W_BIAS) = Q35;
        Q(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.W_BIAS) = Q15;
        Q(obj.INDS_STATE.VEL,obj.INDS_STATE.ACC_BIAS) = Q24;
        Q(obj.INDS_STATE.VEL,obj.INDS_STATE.W_BIAS) = Q25;
        Q(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.ACC_BIAS) = Q44;
        Q(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.W_BIAS) = Q55;
        
        Q(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.VEL) = Q21';
        Q(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.POS) = Q31';
        Q(obj.INDS_STATE.VEL,obj.INDS_STATE.POS) = Q32';
        Q(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.POS) = Q34';
        Q(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.POS) = Q35';
        Q(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.ATTITUDE) = Q15';
        Q(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.VEL) = Q24';
        Q(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.VEL) = Q25';
        
    end
    
    if PARAMS.states.RX_DCB
        Q(obj.INDS_STATE.RX_DCB.INDS,obj.INDS_STATE.RX_DCB.INDS)     =  ...
            eye(length(obj.INDS_STATE.RX_DCB.INDS)) * PARAMS.Q.RX_DCB.^2 * dtKf;
    end
    
    if PARAMS.states.trop
        Q(obj.INDS_STATE.TROP,obj.INDS_STATE.TROP)     =  ...
            eye(length(obj.INDS_STATE.TROP)) * PARAMS.Q.TROP.^2 * dtKf;
    end
    
    if PARAMS.states.iono && strcmp(PARAMS.states.ionoMode,'L1DELAYSTATE')
        indsTec = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
        Q(indsTec,indsTec) = eye(length(indsTec))*PARAMS.Q.L1_IONO.^2*dtKf;
    end
    
    if PARAMS.states.iono && strcmp(PARAMS.states.ionoMode,'TECSTATE')
        indsTec = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
        Q(indsTec,indsTec) = eye(length(indsTec))*PARAMS.Q.L1_IONO.^2*dtKf;
    end
    
    if PARAMS.states.RX_DCB_GLO
        inds_RX_DCB_GLO = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 3 & obj.INDS_STATE.FLEX_STATES_INFO(:,2) == 2);
        Q(inds_RX_DCB_GLO,inds_RX_DCB_GLO) = eye(length(inds_RX_DCB_GLO))*PARAMS.Q.RX_DCB_GLO.^2*dtKf;
    end
    
    if PARAMS.states.RX_DCB_GPS
        inds_RX_DCB_GPS = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 3 & ...
            obj.INDS_STATE.FLEX_STATES_INFO(:,2) == 1);
        Q(inds_RX_DCB_GPS,inds_RX_DCB_GPS) = eye(length(inds_RX_DCB_GPS))*PARAMS.Q.RX_DCB_GPS.^2*dtKf;
    end
    
    if PARAMS.states.MP_CODE
        inds_MP_CODE = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 4);
        Phi_matrix(inds_MP_CODE,inds_MP_CODE) = eye(length(inds_MP_CODE))*exp(-dtKf/PARAMS.other.TAU_MP_CODE);
        Q(inds_MP_CODE,inds_MP_CODE) = eye(length(inds_MP_CODE))*PARAMS.Q.MP_CODE.^2*dtKf;
    end
    
    
    if PARAMS.states.MP_CARR
        inds_MP_CARR = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 5);
        Phi_matrix(inds_MP_CARR,inds_MP_CARR) = eye(length(inds_MP_CARR))*exp(-dtKf/PARAMS.other.TAU_MP_CODE);
        Q(inds_MP_CARR,inds_MP_CARR) = eye(length(inds_MP_CARR))*PARAMS.Q.MP_CARR.^2*dtKf;
    end
    
    indsAmbs = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 1);
    Q(indsAmbs,indsAmbs) = eye(length(indsAmbs))*PARAMS.Q.AMB.^2*dtKf;
    
    % 3. Propagate state estimates using (3.14) noting that only the clock
    % states are non-zero due to closed-loop correction.
    x_est_propagated = zeros(size(obj.state));
    x_est_propagated(obj.INDS_STATE.CLOCK_BIAS,1)   = obj.clockBias+dtKf*obj.clockDrift;
    x_est_propagated(obj.INDS_STATE.CLOCK_DRIFT,1)  = obj.clockDrift;
    x_est_propagated(obj.INDS_STATE.FLEX_STATES,1)  = obj.state(obj.INDS_STATE.FLEX_STATES);
    x_est_propagated(obj.INDS_STATE.TROP,1)         = obj.state(obj.INDS_STATE.TROP);
    x_est_propagated(obj.INDS_STATE.RX_DCB.INDS,1)  = obj.state(obj.INDS_STATE.RX_DCB.INDS);
else
    % inertial updates have not occurred - this is a normal KF time and
    % measurement update
    Phi_matrix = eye(size(cov));
    % velocity -> position
    Phi_matrix(obj.INDS_STATE.POS,obj.INDS_STATE.VEL) = eye(3) * dtKf;
    % clock drift -> clock bias
    Phi_matrix(obj.INDS_STATE.CLOCK_BIAS,obj.INDS_STATE.CLOCK_DRIFT) = dtKf;
    
    
    % 2. Determine approximate system noise covariance matrix using (14.82)
    Q              = zeros(size(cov));
    Q(obj.INDS_STATE.ATTITUDE,obj.INDS_STATE.ATTITUDE) = ...
        eye(3) * PARAMS.Q.gyro_noise_PSD.^2 * dtKf;
    Q(obj.INDS_STATE.POS,obj.INDS_STATE.POS)     =  ...
        eye(3) * PARAMS.Q.POS.^2 * dtKf;
    Q(obj.INDS_STATE.VEL,obj.INDS_STATE.VEL)     =  ...
        eye(3) * PARAMS.Q.VEL.^2 * dtKf;
    Q(obj.INDS_STATE.ACC_BIAS,obj.INDS_STATE.ACC_BIAS) = ...
        eye(3) * PARAMS.Q.ACC_BIAS.^2 * dtKf;
    Q(obj.INDS_STATE.W_BIAS,obj.INDS_STATE.W_BIAS) = ...
        eye(3) * PARAMS.Q.W_BIAS.^2 * dtKf;
    Q(obj.INDS_STATE.CLOCK_BIAS,obj.INDS_STATE.CLOCK_BIAS) = ...
        eye(length(obj.INDS_STATE.CLOCK_BIAS)).*(PARAMS.Q.RXB).^2 * dtKf;
    Q(obj.INDS_STATE.CLOCK_DRIFT,obj.INDS_STATE.CLOCK_DRIFT) = ...
        eye(length(obj.INDS_STATE.CLOCK_BIAS)).*(PARAMS.Q.RXB).^2 * dtKf;
    
    if PARAMS.states.RX_DCB
        Q(obj.INDS_STATE.RX_DCB.INDS,obj.INDS_STATE.RX_DCB.INDS)     =  ...
            eye(length(obj.INDS_STATE.RX_DCB.INDS)) * PARAMS.Q.RX_DCB.^2 * dtKf;
    end
    
    if PARAMS.states.trop
        Q(obj.INDS_STATE.TROP,obj.INDS_STATE.TROP)     =  ...
            eye(length(obj.INDS_STATE.TROP)) * PARAMS.Q.TROP.^2 * dtKf;
    end
    
    if PARAMS.states.iono && strcmp(PARAMS.states.ionoMode,'L1DELAYSTATE')
        indsTec = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
        Q(indsTec,indsTec) = eye(length(indsTec))*PARAMS.Q.L1_IONO.^2*dtKf;
    end
    
    if PARAMS.states.iono && strcmp(PARAMS.states.ionoMode,'TECSTATE')
        indsTec = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 2);
        Q(indsTec,indsTec) = eye(length(indsTec))*PARAMS.Q.L1_IONO.^2*dtKf;
    end
    
    if PARAMS.states.RX_DCB_GLO
        inds_RX_DCB_GLO = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 3 & obj.INDS_STATE.FLEX_STATES_INFO(:,2) == 2);
        Q(inds_RX_DCB_GLO,inds_RX_DCB_GLO) = eye(length(inds_RX_DCB_GLO))*PARAMS.Q.RX_DCB_GLO.^2*dtKf;
    end
    
    if PARAMS.states.RX_DCB_GPS
        inds_RX_DCB_GPS = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 3 & ...
            obj.INDS_STATE.FLEX_STATES_INFO(:,2) == 1);
        Q(inds_RX_DCB_GPS,inds_RX_DCB_GPS) = eye(length(inds_RX_DCB_GPS))*PARAMS.Q.RX_DCB_GPS.^2*dtKf;
    end
    
    if PARAMS.states.MP_CODE
        inds_MP_CODE = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 4);
        
        Phi_matrix(inds_MP_CODE,inds_MP_CODE) = eye(length(inds_MP_CODE))*exp(-dtKf/PARAMS.other.TAU_MP_CODE);
        
        Q(inds_MP_CODE,inds_MP_CODE) = eye(length(inds_MP_CODE))*PARAMS.Q.MP_CODE.^2*dtKf;
    end
    
    if PARAMS.states.MP_CARR
        inds_MP_CARR = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 5);
        Phi_matrix(inds_MP_CARR,inds_MP_CARR) = eye(length(inds_MP_CARR))*exp(-dtKf/PARAMS.other.TAU_MP_CODE);
        Q(inds_MP_CARR,inds_MP_CARR) = eye(length(inds_MP_CARR))*PARAMS.Q.MP_CARR.^2*dtKf;
    end
    
    indsAmbs = obj.INDS_STATE.FLEX_STATES(obj.INDS_STATE.FLEX_STATES_INFO(:,3) == 1);
    Q(indsAmbs,indsAmbs) = eye(length(indsAmbs))*PARAMS.Q.AMB.^2*dtKf;
    
    % Propagate position
    pos = pos+dtKf*vel;
    
    % 3. Propagate state estimates using (3.14) noting that only the clock
    % states are non-zero due to closed-loop correction.
    x_est_propagated = zeros(size(obj.state));
    x_est_propagated(obj.INDS_STATE.CLOCK_BIAS,1)   = obj.clockBias+dtKf*obj.clockDrift;
    x_est_propagated(obj.INDS_STATE.CLOCK_DRIFT,1)  = obj.clockDrift;
    x_est_propagated(obj.INDS_STATE.TROP,1)         = obj.state(obj.INDS_STATE.TROP);
    x_est_propagated(obj.INDS_STATE.FLEX_STATES,1)  = obj.state(obj.INDS_STATE.FLEX_STATES);
    x_est_propagated(obj.INDS_STATE.RX_DCB.INDS,1)  = obj.state(obj.INDS_STATE.RX_DCB.INDS);
    
end

% 4. Propagate state estimation error covariance matrix using (3.46)
cov_propagated = Phi_matrix * (cov + 0.5 * Q) *...
    Phi_matrix' + 0.5 * Q;

% position that is carried around is IMU location- need to adjust for arm
% offset. Also adjust for solid tides
if tSinceInertialUpdate < dtKf
    posRx = pos+R_b_e*PARAMS.IMU_ARM;
    velRx = vel+R_b_e*crossProdMatrix(gyroMeas)*PARAMS.IMU_ARM - ...
        Omega_ie*R_b_e*PARAMS.IMU_ARM*0;
else
    posRx = pos;
    velRx = vel;
end

%% Pull off measurements
measMat = [];

% Wipe off measurements based on the desried measurement masking
obs = navsu.ppp.measMask(obs,PARAMS.measMask);

% Just using all PR measurements
indsMeasPr = find(obs.range.obs ~= 0  & obs.range.ind == 1 );
% indsMeasPr = [];
prnObsMat      = repmat(obs.PRN,size(obs.range.obs,1),1);
constIndObsMat = repmat(obs.constInds,size(obs.range.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnObsMat(indsMeasPr) constIndObsMat(indsMeasPr) ...
    obs.range.sig(indsMeasPr) obs.range.freqs(indsMeasPr) obs.range.obs(indsMeasPr) ...
    1*ones(size(indsMeasPr))]];

% Add carrier phase measurements
indsMeasPh = find(obs.range.obs ~= 0  & obs.range.ind == 2 );
% indsMeasPh = [];
prnObsMat      = repmat(obs.PRN,size(obs.range.obs,1),1);
constIndObsMat = repmat(obs.constInds,size(obs.range.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnObsMat(indsMeasPh) constIndObsMat(indsMeasPh) ...
    obs.range.sig(indsMeasPh) obs.range.freqs(indsMeasPh) obs.range.obs(indsMeasPh) ...
    2*ones(size(indsMeasPh))]];

% Add doppler measurements
indsMeasDop = find(obs.doppler.obs ~= 0);
% indsMeasDop = [];

prnDopMat      = repmat(obs.PRN,size(obs.doppler.obs,1),1);
constIndDopMat = repmat(obs.constInds,size(obs.doppler.obs,1),1);

% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
measMat = [measMat; [prnDopMat(indsMeasDop) constIndDopMat(indsMeasDop) ...
    obs.doppler.sig(indsMeasDop) obs.doppler.freqs(indsMeasDop) ...
    -obs.doppler.obs(indsMeasDop) 3*ones(size(indsMeasDop))]];

nMeas = size(measMat,1);

nSats = length(unique(measMat(:,1:2),'rows'));

if nMeas > 0
    %% Propagate orbit and clock for all measurements
    prnConstInds = sortrows(unique(measMat(:,1:2),'rows'),2);
    
    % transmission time for all satellites
    % rough estimate of travel time- just use geometric range
    [~,ib] = ismember(prnConstInds(:,2),obj.INDS_STATE.CLOCK_BIAS_CONSTS);
    tRx = epoch-obj.clockBias(ib)./c;
    [svPos,svVel] = corrData.propagate(prnConstInds(:,1),prnConstInds(:,2),tRx);
    
    % Might be missing some precise data- check and remove if so
    if any(isnan(svPos(:,1)))
        indsNan = find(isnan(svPos(:,1)));
        
        prnConstIndsNan = prnConstInds(indsNan,:);
        
        % Remove the associated measurements
        indsMeasRemove = find(ismember(measMat(:,1:2),prnConstIndsNan,'rows'));
        measMat(indsMeasRemove,:) = [];
        
        prnConstInds(indsNan,:) = [];
        [~,ib] = ismember(prnConstInds(:,2),obj.INDS_STATE.CLOCK_BIAS_CONSTS);
        tRx = epoch-obj.clockBias(ib)./c;
        
        svPos(indsNan,:) = [];
        
        nMeas = size(measMat,1);
    end
    
    satBias = -c*corrData.clock(prnConstInds(:,1),prnConstInds(:,2),tRx);
    gRangeSv = sqrt(sum((obj.pos'-svPos).^2,2));
    
    % Need rough estimate of the receiver clock bias in case of reset
    measMatPr = measMat(measMat(:,6) == 1,:);
    % Pull one pseudorange for each satellite
    [~,losIndPr] = ismember(measMatPr(:,1:2),prnConstInds,'rows');
    
    % bRxi = nanmedian(measi(measInfoi(:,4) == 1,1)- sqrt(sum((svPosRot-repmat(usrPos',size(svPosRot,1),1)).^2,2))-bSati(sIndsMap(measInfoi(:,4) == 1)));
    if isempty(measMatPr) || 1
        bRxi = obj.clockBias;
        bRxi = x_est_propagated(obj.INDS_STATE.CLOCK_BIAS,1);
    else
        bRxi = nanmedian(measMatPr(:,5)-gRangeSv(losIndPr)-satBias(losIndPr));
        obj.clockBias(:) = bRxi;
    end
    x_est_propagated(obj.INDS_STATE.CLOCK_BIAS,1)   = bRxi;
    
    epoch = epoch-obj.clockBias(ib)./c;
    
    satBias = -c*corrData.clock(prnConstInds(:,1),prnConstInds(:,2),epoch);
    [svPos,svVel] = corrData.propagate(prnConstInds(:,1),prnConstInds(:,2),epoch);
    gRangeSv = sqrt(sum((obj.pos'-svPos).^2,2));
    
    tTx = epoch-gRangeSv./c;
    
    [svPos,svVel] = corrData.propagate(prnConstInds(:,1),prnConstInds(:,2),tTx);
    travelTime = epoch-tTx;
    svPosRot = navsu.ppp.models.earthRotTravelCorr(travelTime,svPos);
    
    rxDrift   = obj.clockDrift(ib);
    rxBias    = obj.clockBias(ib);
    
    % elevation and azimuth for each LOS
    [el,az] = navsu.geo.pos2elaz(obj.pos',svPosRot);
    
    %% Various range effects
    % tropo delay for each LOS
    % trop = saastamoinen_model_SU(llhi(1), llhi(2), llhi(3), el*180/pi);
    % [trop0,m,tropDataSave] = tropo_error_correction_unb3(el,h,lat,doy);
    doy = navsu.time.jd2doy(navsu.time.epochs2jd(epoch));
    [trop,m,~] = navsu.ppp.models.tropDelay(el*180/pi,az*180/pi, llhi(:,3), llhi(:,1), llhi(:,2), doy, PARAMS, [],[],epoch);
    
    % TEC for each LOS
    if any(measMat(:,3) < 100 & measMat(:,6) < 3) %&& strcmp(PARAMS.states.ionoMode,'TEC')
        [~,~,tecSlant] = corrData.ionoDelay(epoch,llhi,'az',az,'el',el);
    else
        tecSlant = zeros(size(prnConstInds,1),1);
    end
    
    % solid tide adjustment to position
    [~,stRangeOffset] = navsu.ppp.models.solidTide(epoch(1),pos,'svPos',svPosRot);
    
    % relativistic corrections
    relClockCorr = 2/c^2.*sum(svPos.*svVel,2)*c;
    relRangeCorr = navsu.ppp.models.relRangeCorr(svPos',pos',PARAMS);
    
    % Carrier phase windup
    [~,ib] = ismember(prnConstInds,obj.phWind.PrnConstInd,'rows');
    phWind = navsu.ppp.models.carrierPhaseWindupGGM(epoch(1), repmat(pos',size(svPosRot,1)), svPosRot, obj.phWind.phaseOffset(ib));
    obj.phWind.phaseOffset(ib) = phWind; % need to update the phase windup object
    
    
    
    % Geometry matrix
    A = (svPosRot-pos')./sqrt(sum((pos'-svPosRot).^2,2));
    
    % geometric range
    gRange = sqrt(sum((svPosRot-posRx').^2,2));
    dVel   = dot([+svVel-velRx']',-A')';
    dVel2 = velRx'-svVel;
    
    
    %%
    [~,losInds] = ismember(measMat(:,[1 2]),prnConstInds,'rows');
    [~,indAmbStates] = ismember([measMat(:,[1 2 3]) ones(size(measMat,1),1)],obj.INDS_STATE.FLEX_STATES_INFO(:,[1 2 4 3]),'rows');
    % indAmbStates = obj.INDS_STATES.FLEX_STATES(ia);
    % indIono = obj.INDS_STATE.FLEX_STATES(ismember(obj.INDS_STATE.FLEX_STATES_INFO(:,1:3),[prni constIndi 2],'rows'));
    [~,indIonos] =ismember([measMat(:,1:2) 2*ones(size(measMat,1),1)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:3),'rows');
    [~,indGloDcbs] =ismember([measMat(:,1:2) 3*ones(size(measMat,1),1) measMat(:,3)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    % [~,indGpsDcbs] =ismember([measMat(:,1:2) 3*ones(size(measMat,1),1) measMat(:,3)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    
    [~,indMpCodes] =ismember([measMat(:,1:2) 4*ones(size(measMat,1),1) measMat(:,3)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    [~,indMpCarrs] =ismember([measMat(:,1:2) 5*ones(size(measMat,1),1) measMat(:,3)],obj.INDS_STATE.FLEX_STATES_INFO(:,1:4),'rows');
    
    % build each measurement
    pred_meas = zeros(nMeas,1);
    H         = zeros(nMeas,nState);
    R         = zeros(nMeas,1);
    for idx = 1:nMeas
        measTypei = measMat(idx,6);% 1 PRN | 2 const | 3 signal (sf/df) | 4 freq | 5 meas | 6 (1=pr,2=ph,3=dop)
        prni      = measMat(idx,1);
        constIndi = measMat(idx,2);
        sigi      = measMat(idx,3);
        freqi     = measMat(idx,4);
        measi     = measMat(idx,5);
        
        losInd = losInds(idx);
        weighti = 1 ./ (sin(el(losInd)).^2);
        
        if nSats < 0
            %             weighti = weighti * 10000;
            sigAddRange = 100^2;
            sigAddDoppler = 100^2;
        else
            sigAddRange = 0;
            sigAddDoppler = 0;
        end
        
        switch measTypei
            case 1
                %% code
                if sigi < 100 && strcmp(PARAMS.states.ionoMode,'TEC')
                    % If it's single frequency, need to include iono correction
                    ionoCorri = -tecSlant(losInd)*40.3*10^15./freqi.^2;
                    
                elseif sigi < 100 && strcmp(PARAMS.states.ionoMode,'L1DELAYSTATE')
                    
                    ionoCorrModel = -tecSlant(losInd)*40.3*10^15./freqi.^2;
                    
                    indIono = obj.INDS_STATE.FLEX_STATES(indIonos(idx));
                    delayL1i = x_est_propagated(indIono);
                    
                    hi = -(1575.42e6).^2./freqi.^2;
                    H(idx,indIono) = hi;
                    ionoCorri = delayL1i*hi+ionoCorrModel;
                else
                    ionoCorri = 0;
                end
                
                if PARAMS.states.RX_DCB && ~(PARAMS.states.RX_DCB_GLO && constIndi == 2) ...
                        && ~(PARAMS.states.RX_DCB_GPS && constIndi == 1)
                    indRxDcb = find(obj.INDS_STATE.RX_DCB.sig == sigi & obj.INDS_STATE.RX_DCB.constInds == constIndi, 1);
                    if ~isempty(indRxDcb)
                        rxDcb = x_est_propagated(obj.INDS_STATE.RX_DCB.INDS(indRxDcb));
                        H(idx,obj.INDS_STATE.RX_DCB.INDS(indRxDcb)) = 1;
                    else
                        rxDcb = 0;
                    end
                else
                    rxDcb = 0;
                end
                
                if PARAMS.states.RX_DCB_GPS && constIndi == 1 && indGloDcbs(idx)~= 0
                    indDcbGlo = obj.INDS_STATE.FLEX_STATES(indGloDcbs(idx));
                    dcbGpsi = x_est_propagated(indDcbGlo);
                    H(idx,indDcbGlo) = 1;
                else
                    dcbGpsi = 0;
                end
                
                if PARAMS.states.RX_DCB_GLO && constIndi == 2 && indGloDcbs(idx)~= 0
                    indDcbGlo = obj.INDS_STATE.FLEX_STATES(indGloDcbs(idx));
                    dcbGloi = x_est_propagated(indDcbGlo);
                    H(idx,indDcbGlo) = 1;
                else
                    dcbGloi = 0;
                end
                
                if PARAMS.states.MP_CODE
                    indMpCode = obj.INDS_STATE.FLEX_STATES(indMpCodes(idx));
                    mpCodei   = x_est_propagated(indMpCode);
                    H(idx,indMpCode) = 1;
                else
                    mpCodei = 0;
                end
                
                if PARAMS.states.trop
                    dtrop = m(losInd)*x_est_propagated(obj.INDS_STATE.TROP);
                    H(idx,obj.INDS_STATE.TROP) = m(losInd);
                else
                    dtrop = 0;
                end
                
                pred_meas(idx) = gRange(losInd)+satBias(losInd)+rxBias(losInd)+trop(losInd)+dtrop+...
                    stRangeOffset(losInd)+relClockCorr(losInd)+relRangeCorr(losInd)+ionoCorri+rxDcb+...
                    dcbGloi+dcbGpsi+mpCodei;
                
                H(idx,obj.INDS_STATE.POS)        = A(losInd,:);
                H(idx,obj.INDS_STATE.CLOCK_BIAS(constIndi)) = 1;
                
                if sigi < 100
                    R(idx,idx) = weighti.*3^2+sigAddRange;
                else
                    R(idx,idx) = weighti.*3^2+sigAddRange;
                end
                
            case 2
                %% carrier
                if sigi < 100 && strcmp(PARAMS.states.ionoMode,'TEC')
                    % If it's single frequency, need to include iono correction
                    ionoCorri = +tecSlant(losInd)*40.3*10^15./freqi.^2;
                    
                elseif sigi < 100 && strcmp(PARAMS.states.ionoMode,'L1DELAYSTATE')
                    ionoCorrModel = -tecSlant(losInd)*40.3*10^15./freqi.^2;
                    
                    indIono = obj.INDS_STATE.FLEX_STATES(indIonos(idx));
                    delayL1i = x_est_propagated(indIono);
                    
                    hi = (1575.42e6).^2./freqi.^2;
                    H(idx,indIono) = hi;
                    
                    ionoCorri = hi*delayL1i+ionoCorrModel;
                else
                    ionoCorri = 0;
                end
                
                if PARAMS.states.MP_CARR
                    indMpCarr = obj.INDS_STATE.FLEX_STATES(indMpCarrs(idx));
                    mpCarri   = x_est_propagated(indMpCarr);
                    H(idx,indMpCarr) = 1;
                else
                    mpCarri = 0;
                end
                
                if PARAMS.states.trop
                    dtrop = m(losInd)*x_est_propagated(obj.INDS_STATE.TROP);
                    H(idx,obj.INDS_STATE.TROP) = m(losInd);
                else
                    dtrop = 0;
                end
                
                indAmbState = obj.INDS_STATE.FLEX_STATES(indAmbStates(idx));
                
                ambEst = x_est_propagated(indAmbState);
                
                % Carrier phase windup
                phWindi = phWind(losInd)*c/freqi;
                
                pred_meas(idx) = gRange(losInd)+satBias(losInd)+rxBias(losInd)+...
                    trop(losInd)+dtrop+stRangeOffset(losInd)+relClockCorr(losInd)+...
                    relRangeCorr(losInd)+ionoCorri+ambEst+phWindi+mpCarri;
                
                H(idx,obj.INDS_STATE.POS)        = A(losInd,:);
                H(idx,obj.INDS_STATE.CLOCK_BIAS(constIndi)) = 1;
                H(idx,indAmbState)              = 1;
                
                if sigi < 100
                    R(idx,idx) = weighti.*0.03^2+sigAddRange;
                else
                    R(idx,idx) = weighti.*0.03^2+sigAddRange;
                end
                
            case 3
                % doppler
                if 0
                    [hVel,hAtt,hBias] = dopplerSensitivity(vel,svVel(losInd,:)',rxDrift(losInd),A(losInd,:),...
                        gyroMeas,R_b_e,PARAMS.IMU_ARM,Omega_ie);
                    
                    %                 pred_meas(idx) = dopplerModel(vel,svVel(losInd,:)',rxDrift(losInd),A(losInd,:),...
                    %                     gyroMeas,R_b_e,PARAMS.IMU_ARM,Omega_ie);
                    pred_meas(idx) = -dot(dVel2(losInd,:),-A(losInd,:))-rxDrift(losInd);
                    
                    H(idx,obj.INDS_STATE.VEL)                    = -hVel;
                    H(idx,obj.INDS_STATE.ATTITUDE)               = -hAtt;
                    H(idx,obj.INDS_STATE.CLOCK_DRIFT(constIndi)) = -1;
                    H(idx,obj.INDS_STATE.W_BIAS)                 = hBias;
                else
                    
                    pred_meas(idx) = -dot(dVel2(losInd,:),-A(losInd,:))-rxDrift(losInd);
                    
                    H(idx,obj.INDS_STATE.VEL) = -A(losInd,:);
                    H(idx,obj.INDS_STATE.CLOCK_DRIFT(constIndi)) = -1;
                end
                
                R(idx,idx) = weighti*0.05^2+sigAddDoppler;
                
                measMat(idx,5) = -measMat(idx,5);
        end
    end
    
    % Remove measurements from satellites that are too low
    indsElLow = find(el(losInds)*180/pi < PARAMS.elMask);
    if ~isempty(indsElLow)
        measMatRemovedLow = measMat(indsElLow,:);
        
        H(indsElLow,:) = [];
        R(indsElLow,:) = [];
        R(:,indsElLow) = [];
        measMat(indsElLow,:) = [];
        pred_meas(indsElLow) = [];
    else
        measMatRemovedLow = zeros(0,size(measMat,2));
    end
else
    % No GNSS measurements
    % build each measurement
    pred_meas = zeros(0,1);
    H         = zeros(0,nState);
    R         = zeros(0,1);
    
end

if nMeas  == 0
    el = [];
    az = [];
    prnConstInds = zeros(0,2);
    measMatRemovedLow = zeros(0,2);
end

%% Pseudomeasurements
if PARAMS.measUse.noVertVel
    % No vertical velocity constraint (vehicle only moves forward)
    pseudoMeasi = 0;
    vertVeli = [0 0 1]*R_b_e'*velRx;
    
    ri = 1^2;
    Hi = zeros(1,size(H,2));
    Hi(1,obj.INDS_STATE.VEL) = -[0 0 1]*R_b_e';
    
    measMati = [0 0 0 0 pseudoMeasi 4];
    
    R(size(R,1)+1,size(R,1)+1) = ri;
    H = [H; Hi];
    measMat = [measMat; measMati];
    pred_meas = [pred_meas; vertVeli];
    
    pseudoMeasi = 0;
    vertVeli = [0 1 0]*R_b_e'*velRx;
    
    ri = 1^2;
    Hi = zeros(1,size(H,2));
    Hi(1,obj.INDS_STATE.VEL) = -[0 1 0]*R_b_e';
    
    measMati = [0 0 0 0 pseudoMeasi 4];
    
    R(size(R,1)+1,size(R,1)+1) = ri;
    H = [H; Hi];
    measMat = [measMat; measMati];
    pred_meas = [pred_meas; vertVeli];
    
end


if ~isempty(pos0) && 0
    [dposEnu,Renui] = XYZ2ENU([pos-pos0]',llhi(1)*pi/180,llhi(2)*pi/180);
    
    ri = 0.01^2;
    Hi = zeros(1,size(H,2));
    Hi(1,obj.INDS_STATE.POS) = Renui(3,:);
    pseudoMeasi = 0;
    measMati = [0 0 0 0 pseudoMeasi 4];
    
    dposUpi     = -dposEnu(3);
    
    R(size(R,1)+1,size(R,1)+1) = ri;
    H = [H; Hi];
    measMat = [measMat; measMati];
    pred_meas = [pred_meas; dposUpi];
end

nMeas = size(H,1);
%%
if nMeas > 0
    % Measurement were available- do the update.
    [H,delta_z,residsPost,K,measMat,~,measMatRemoved,R] = ...
        navsu.ppp.measUpdateExclude(H,cov_propagated,R,measMat,pred_meas,PARAMS);
    
    % 9. Update state estimates
    x_est_new = x_est_propagated + K * delta_z;
    
    % Update covariance
    cov = (eye(nState) - K * H) * cov_propagated;
    
else
    % No measurement update
    x_est_new = x_est_propagated;
    
    cov = cov_propagated;
    
    measMatRemoved = zeros(0,6);
    residsPost = [];
end

% CLOSED-LOOP CORRECTION

% Correct attitude, velocity, and position using (14.7-9)
R_b_e = (eye(3) - navsu.geo.crossProdMatrix(x_est_new(obj.INDS_STATE.ATTITUDE))) * R_b_e;

if tSinceInertialUpdate < dtKf
    velRx = velRx - x_est_new(obj.INDS_STATE.VEL);
    posRx = posRx - x_est_new(obj.INDS_STATE.POS);
    
    pos = posRx-R_b_e*PARAMS.IMU_ARM;
    vel = velRx-R_b_e*crossProdMatrix(gyroMeas)*PARAMS.IMU_ARM + ...
        Omega_ie*R_b_e*PARAMS.IMU_ARM;
else
    vel = vel - x_est_new(obj.INDS_STATE.VEL);
    pos = pos - x_est_new(obj.INDS_STATE.POS);
    
end

obj.clockBias  = x_est_new(obj.INDS_STATE.CLOCK_BIAS);
obj.clockDrift = x_est_new(obj.INDS_STATE.CLOCK_DRIFT);


% Update IMU bias and GNSS receiver clock estimates
imuBiasStates = imuBiasStates + x_est_new([obj.INDS_STATE.ACC_BIAS obj.INDS_STATE.W_BIAS]);

% put updated values into object
obj.R_b_e = R_b_e;
obj.vel   = vel;
obj.pos   = pos;
obj.imuBiasStates = imuBiasStates;
obj.cov  = cov;
obj.state = x_est_new;
obj.posPrevTc = pos;

% Deal with resets if any of the removed measurements were carrier phases
if ~isempty(measMatRemoved)
    for jdx = 1:size(measMatRemoved,1)
        if measMatRemoved(jdx,6) == 2 %|| measMatRemoved(jdx,6) == 1
            % carrier phase- reset ambiguity by just removing the state
            obj.removeFlexState([measMatRemoved(jdx,[1 2]) 1 measMatRemoved(jdx,3)] );
        end
    end
end

%% Update the list of measurements used (if any haven't been included yet)
satsUsedi = unique(measMat(:,1:2),'rows');
obj.allSatsSeen = sortrows(unique([measMat(:,1:2); obj.allSatsSeen],'rows'),2);

%% Save for output
if ~isempty(outStruc)
    [rangeResids,doppResids,elFull,azFull] = outStruc.saveResids(measMat,residsPost,epoch0,el,az,prnConstInds);
    
    outStruc.saveMeasRemoved(epoch0,measMatRemovedLow,measMatRemoved);
end

obj.resids.epoch = epoch0;
obj.resids.range = rangeResids;
obj.resids.doppler = doppResids;
obj.resids.el      = elFull;
obj.resids.az      = azFull;

end






















