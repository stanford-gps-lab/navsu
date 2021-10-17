classdef CarrierSmoother < matlab.mixin.Copyable
    %CarrierSmoother Class to perform carrier smoothing.
    %   Performs smoothing recursively leveraging previous signals. One
    %   object should be instantiated per frequency.
    
    properties
        lastCarrierPhase    % Value of last carrier phase measurement in m
        tLastSmoothing      % Time of last carrier phase measurement in epoch
        lastSmoothedPr      % Value of last smoothed pseudorange in m
        nSmoothingEpochs    % Number of epochs that each range has been smoothened
        lastPrVar           % variance of last smoothened pseudorange value
        tMax                % maximum time period of smoothing to avoid 
                            % code-carrier divergence
    end
    
    methods
        function obj = CarrierSmoother(numSat, tMax)
            %CarrierSmoother(numSat, frequency) Create a smoothing object.
            %   Instantiate a carrier smoothing object for a given number
            %   of satellites.
            
            if nargin < 2
                tMax = 100;
            end
            if nargin < 1
                numSat = 0;
            end
            
            obj.lastCarrierPhase    = NaN(numSat, 1);
            obj.tLastSmoothing      = NaN(numSat, 1);
            obj.lastSmoothedPr      = NaN(numSat, 1);
            obj.nSmoothingEpochs    = zeros(numSat, 1);
            obj.lastPrVar           = zeros(numSat, 1);
            
            obj.tMax = tMax;
            
        end
        
        function [prSmoothed, prVarRcvr] = smoothen(obj, ...
                satIds, epoch, code, carrier, lockTime, codeVar, carrierVar)
            % Perform carrier smoothing for given code, carrier,
            % observations with a given carrier lock time at a given epoch.
            
            % set default inputs
            if nargin < 7 || isempty(codeVar)
                codeVar = NaN(size(code));
            end
            if nargin < 8 || isempty(carrierVar)
                carrierVar = NaN(size(carrier));
            end
            
            % get difference in carrier phase by which to update
            carrierDelta = carrier - obj.lastCarrierPhase(satIds);
            
            
            % first step: get smoothing constant M for each satellite
            M = obj.nSmoothingEpochs(satIds); % how many epochs have already been smoothed
            
            % increase by 1 for every sat with finite carrier phase and
            % maintained carrier lock, else set to 1
            dt = epoch - obj.tLastSmoothing(satIds);
            M = 1 + M .* (isfinite(carrierDelta) ...
                        & lockTime > dt ...
                        & lockTime < obj.tMax);
            M = min(M, obj.tMax./dt);
            % -> M is now the number of epochs across which to smoothen
            
            
            % second step: actually do the smoothing
            prUpdate = (M-1)./M .* (obj.lastSmoothedPr(satIds) + carrierDelta);
            prUpdate(M == 1) = 0;
            prSmoothed = code ./ M + prUpdate;
            
            % update sigma accordingly
            varUpdate = carrierVar + obj.lastPrVar(satIds);
            varUpdate(M == 1) = 0;
            prVarRcvr = (codeVar + (M-1).^2 .* varUpdate) ./ M.^2;
            
            
            % third step: update stored values
            obj.lastCarrierPhase(satIds)    = carrier;
            obj.lastSmoothedPr(satIds)      = prSmoothed;
            obj.tLastSmoothing(satIds)      = epoch;
            obj.nSmoothingEpochs(satIds)    = M;
            obj.lastPrVar(satIds)           = prVarRcvr;
        end
        
        
    end
end

