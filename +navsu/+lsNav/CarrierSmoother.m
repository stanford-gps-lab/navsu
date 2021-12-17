classdef CarrierSmoother < matlab.mixin.Copyable
    %CarrierSmoother Class to perform carrier smoothing.
    %   Performs smoothing recursively leveraging previous signals. One
    %   object should be instantiated per frequency.
    %   
    %   Properties:
    %   signal              Rinex identifier of the smoothed signals
    %   const               constellation indicator
    %   lastCarrierPhase    Value of last carrier phase measurement in m
    %   tLastSmoothing      Time of last carrier phase measurement in epoch
    %   lastSmoothedPr      Value of last smoothed pseudorange in m
    %   nSmoothingEpochs    Number of epochs that each range has been smoothened
    %   lastPrVar           Variance of last smoothened pseudorange value
    %   tMax                maximum time period of smoothing to avoid
    %                       code-carrier divergence
    
    properties
        signal(1, :) char = '1C'
        const(1, 1) double {mustBePositive, mustBeFinite} = 1
        tMax(1, 1) double = 100
        lastCarrierPhase double
        tLastSmoothing double
        lastSmoothedPr double
        nSmoothingEpochs double
        lastPrVar double       
    end
    
    methods
        function obj = CarrierSmoother(signal, const, numSat, tMax)
            %CarrierSmoother(signal, const, numSat, frequency)
            %   Instantiate a carrier smoothing object for a given number
            %   of satellites.
            
            % read and store inputs
            
            if nargin > 0
                obj.signal = signal;
            end
            if nargin > 1
                obj.const = const;
            end
            if nargin < 3
                numSat = 0;
            end
            if nargin > 3
                obj.tMax = tMax;
            end
            
            obj.lastCarrierPhase    = NaN(numSat, 1);
            obj.tLastSmoothing      = NaN(numSat, 1);
            obj.lastSmoothedPr      = NaN(numSat, 1);
            obj.nSmoothingEpochs    = zeros(numSat, 1);
            obj.lastPrVar           = zeros(numSat, 1);
                        
        end
        
        function [prSmoothed, prVarRcvr] = smoothen(obj, ...
                satIds, epoch, obsData, SigRNM)
            % Perform carrier smoothing for given code, carrier,
            % observations with a given carrier lock time at a given epoch.

            % set default inputs
            if nargin < 5 || isempty(SigRNM)
                SigRNM.code = NaN(size(obsData.code));
                SigRNM.carrier = NaN(size(obsData.carrier));
            end
            
            % % could consider doppler as check on cycle slips
            % lCp = [obj.CS(1).lastCarrierPhase(satIds) obj.CS(2).lastCarrierPhase(satIds)];
            % lCpE = [obj.CS(1).tLastSmoothing(satIds) obj.CS(2).tLastSmoothing(satIds)];
            % % get rate of change of carrier phase in cycles/sec
            % carrierDelta = (carrierMeas(:, 1:2) - lCp) ./ (epoch - lCpE) ...
            %               .* freq / navsu.constants.c;
            % carrierVSdoppler = carrierDelta + obsData.doppler(~noPr,  :);
            % % should be close to 0

            % first step: get smoothing constant M for each satellite
            M = obj.getSmoothingConstant(satIds, epoch, obsData.tLock);
            % -> M is now the number of epochs across which to smoothen
            
            % get difference in carrier phase by which to update
            carrierDelta = obsData.carrier - obj.lastCarrierPhase(satIds);
            
            % second step: actually do the smoothing
            prUpdate = (M-1)./M .* (obj.lastSmoothedPr(satIds) + carrierDelta);
            badUpdate = ~isfinite(prUpdate);
            prUpdate(badUpdate) = 0;
            M(badUpdate) = 1;
            prSmoothed = obsData.code ./ M + prUpdate;
            
            % update sigma accordingly
            varUpdate = SigRNM.carrier + obj.lastPrVar(satIds);
            varUpdate(isnan(varUpdate)) = 0;
            prVarRcvr = (SigRNM.code + (M-1).^2 .* varUpdate) ./ M.^2;
            
            
            % third step: update stored values
            obj.lastCarrierPhase(satIds)    = obsData.carrier;
            obj.lastSmoothedPr(satIds)      = prSmoothed;
            obj.tLastSmoothing(satIds)      = epoch;
            obj.nSmoothingEpochs(satIds)    = M;
            obj.lastPrVar(satIds)           = prVarRcvr;
        end
        
        function reset(obj, satIds)
            % Resets the smoother to initial conditions. Can be run for
            % certain satellites only.
            %   
            %   obj.reset;
            %   obj.reset(satIds);

            if nargin < 2
                satIds = true(size(obj.lastCarrierPhase));
            end

            obj.lastCarrierPhase(satIds) = NaN;
            obj.tLastSmoothing(satIds) = NaN;
            obj.lastSmoothedPr(satIds) = NaN;
            obj.nSmoothingEpochs(satIds) = 0;
            obj.lastPrVar(satIds) = 0;

        end
    end

    methods (Access = protected)
        function M = getSmoothingConstant(obj, satIds, epoch, lockTime)
            % Calculate the smoothing constant M for each satellite.
            % M is equal to the number of epochs over which the
            % measurements have been smoothened. Larger M result in
            % stronger weighting of the carrier phase.

            dt = epoch - obj.tLastSmoothing(satIds);
            if any(dt < 0)
                disp(['Resetting ', navsu.svprn.convertConstIndName(obj.const), ...
                      ' carrier smoother for ', obj.signal, ...
                      ' signal due to negative timestep.']);
                obj.reset(satIds);
                % recompute dt
                dt = epoch - obj.tLastSmoothing(satIds);
            end

            % get initial smoothing constant M for each satellite
            M = max(1, obj.nSmoothingEpochs(satIds));
            
            % increase by 1 for every sat with maintained carrier lock,
            % else set to 1
            M = 1 + M .* (lockTime > dt & dt > 0);
            % limit to max smoothing epoch
            M = min(M, obj.tMax./dt);
        end
    end
end

