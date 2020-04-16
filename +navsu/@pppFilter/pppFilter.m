classdef pppFilter < matlab.mixin.Copyable


    properties

        state   % the state of the filter -> should be of type pppanal.ppp.State
        cov     % covariance of the state
        
        pos                % ECEF position
        vel                % ECEF velocity
        R_b_e              % DCM from body to ECEF
        imuBiasStates      % imu bias states
        clockBias          % receiver clock bias(es)
        clockDrift         % receiever clock drift(s)
        carrierAmbiguities % carrier phase ambiguity estimates
        
        INDS_STATE % state indexing
        
        lastAccMeas % most recent acceleration measurements
        lastGyroMeas % most recent gyro measurement
        epochLastInertialUpdate % time of the last inertial update
        epochLastGnssUpdate      % time of the last GNSS update
        
        posPrevTc % position at previous tight coupling update
        
        % contains latest geometry free combinations
        cycleSlipInfo = struct('gFree',zeros(0,1),'epochLastGFree',zeros(0,1),'measInfoGFree',zeros(0,4));
        
        % phase windup
        phWind = struct('phaseOffset',zeros(0,1),'PrnConstInd',zeros(0,2));
                
        StateMap % mapping matrix- indicates what the state is in each position of the state and covariance matrix
        
        % all satellites used in the solution- useful for solution
        % separation :)
        allSatsSeen
        
        initialized = false % whether or not the fitler has been initialized
        
    end


    methods
        function obj = EKF()
            if nargin < 1
                % if no values are provided (this is the only option right
                % now)
                obj.QVals = [];
                obj.RVals = [];
            else
                % if values are provided
                obj.QVals = QVals;
                obj.QVals = QVals;
            end
        end
    end


    % function signatures
    methods
        complete = initialize(obj,varargin)
        inertialNavEquationsEcef(obj,epoch,accMeasi,gyroMeasi);
        lcUpdate(obj,epoch,posMeasi,velMeasi,accMeasi,gyroMeasi,PARAMS)
        [extraInputs,measMatRemoved,measMatLow] = tcUpdate(obj,epoch,gnssMeas,...
            accMeasi,gyroMeasi,PARAMS,outStruc,extraInputs,measMask)
        manageStates(obj,epoch,gnssMeas,PARAMS,outStruc);
        measRemoved = checkCycleSlips(obj,epoch,gnssMeas,PARAMS);        
        removeFlexState(obj,measInfoRemove);
    end


end















