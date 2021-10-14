classdef DFMCnavigationEngine < matlab.mixin.Copyable
    %DFMCnavigationEngine Dual frequency, multi constellation navigation.
    %   Computes navigation solutions from dual frequency, multi
    %   constellation GNSS measurements. Can compute position, velocity and
    %   time solutions.
    
    
    properties
        velocity (3,1)   % current velocity estimate in ECEF in (m/s)
        tBias      % time bias for each constellation in (s)
        tBiasRate   % time bias rate for each constellation in (s/s)
        tPosSol = NaN; % time of last position solution in sec since first GPS epoch
        numConsts   % number of included constellations
        numSats     % number of included satellites
        constellations % char vector indicating the useable constellations
        navCov      % Position, time bias covariance matrix
        rateCov     % Velocity, time rate covariance matrix
        satAcc      % Variance of satellite position
        satPRN      % PRN of each satellite
        satConstId  % constellation id of each satellite
        constLetter % one letter identifier of each constellation: GRECJ
        satTGD      % Timing group delay of each satellite in sec (see IS﻿20.3.3.3.3.2)
        satEl       % elevation of satellites in deg
        satAz       % azimuth of satellites in deg
        freqMap     % lookup table for frequency of each signal
        CS (3, 1) navsu.lsNav.CarrierSmoother % array of objects for carrier 
        % smoothing. One smoother for each freq and one for Dual Freq.
        elevMask = 15*pi/180;    % elevation mask angle
    end
    
    properties (Dependent)
        position    % current position estimate in ECEF coordinates in (m)
        positionLLH % position estimate in lat, lon, height (deg, deg, m)
        losVectors  % line of sight vectors user to satellite
        theoranges  % ranges user to satellite (length of losVectors)
        satPos      % satellite position at time of broadcast
        satVel      % satellite velocity at time of broadcast
        Beph        % broadcast ephemeris, output of navsu.readfiles.loadRinexNav
        Peph        % precise ephemeris object
    end
    
    properties (Access = private, Hidden = true)
        % internal memory
        internal_position (3,1)
        internal_losVectors
        internal_theoranges
        internal_satPos = zeros(0, 3);
        internal_satVel = zeros(0, 3);
        internal_satClkBias
        internal_satClkRate
        internal_satEpoch       % epoch of last orbit propagation
        internal_Beph struct
        internal_Peph navsu.svOrbitClock
    end
    
    methods % SET & GET methods
        
        function sP = get.satPos(obj)
            % retrieve from memory
            sP = obj.internal_satPos;
        end
        function set.satPos(obj, satPositions)
            % set stored value
            obj.internal_satPos = satPositions;
            
            % also recompute the line of sight vectors
            obj.losVectors = obj.satPos - obj.position';
        end
        
        function sV = get.satVel(obj)
            % retrieve from memory
            sV = obj.internal_satVel;
        end
        function set.satVel(obj, satVel)
            obj.internal_satVel = satVel;
        end
        
        function los = get.losVectors(obj)
            % retrieve from memoty
            los = obj.internal_losVectors;
        end
        function set.losVectors(obj, los)
            % set stored value
            obj.internal_losVectors = los;
            % also recompute the theoranges
            obj.theoranges = sqrt(sum(obj.losVectors.^2, 2));
            
            % also recompute az, el
            [obj.satEl, obj.satAz] = navsu.geo.pos2elaz(obj.position', obj.satPos);
        end
        
        function tr = get.theoranges(obj)
            % retrieve from memory
            tr = obj.internal_theoranges;
        end
        function set.theoranges(obj, ranges)
            % set stored value
            obj.internal_theoranges = ranges;
        end
        
        function eph = get.Beph(obj)
            % retrieve from memory
            eph = obj.internal_Beph;
        end
        function set.Beph(obj, eph)
            % set the property
            obj.internal_Beph = eph;
            
            % also initialize internal sat memory properties
            obj.initializeEngine(~structfun(@isempty, eph)');
        end
        
        function eph = get.Peph(obj)
            % retrieve from memory
            eph = obj.internal_Peph;
        end
        function set.Peph(obj, eph)
            % set the property
            % this is a start to implement it
            obj.internal_Peph = eph;
            % also initialize internal sat memory properties
            obj.initializeEngine(eph.settings.constUse);
        end
        
        function llh = get.positionLLH(obj)
            % convert the stored ECEF position to lat, lon, height w.r.t.
            % the WGS84 ellipsoid.
            llh = navsu.geo.xyz2llh(obj.position')';
        end
        
        function pos = get.position(obj)
            
            pos = obj.internal_position;
        end
        function set.position(obj, pos)
            % save value to memory
            obj.internal_position = pos;
            
            % also recompute the line of sight vectors
            obj.losVectors = obj.satPos - obj.position';
        end
        
        
    end
    
    methods
        function obj = DFMCnavigationEngine(eph, x0)
            % DFMCnavigationEngine GNSS navigation engine for DF, MC.
            %   Provides recursive least squares estimation for dual
            %   frequency (DF), multi constellation (MC) GNSS measurements.
            %   
            %   DFMCnavigationEngine(eph, x0)
            %   
            %   Can be initialized without inputs or with up to two inputs.
            %   
            %   Inputs:
            %   eph     Ephemeris information. Can be navsu.svOrbitClock
            %           precise orbit product class or ephemeris struct
            %           that is result of navsu.readfiles.loadRinexNav.
            %   x0      Initial position estimate to accelerate computation
            %           of first fix. Default is [0 0 0]'.
            
            if nargin > 0
                if isa(eph, 'navsu.svOrbitClock')
                    % got precise ephemeris. This one will be preferred.
                    obj.Peph = eph;
                elseif isstruct(eph)
                    % work with broadcast ephemeris
                    obj.Beph = eph;
                end
            end
            
            if nargin > 1
                obj.position = x0;
            else
                obj.position = zeros(3, 1);
            end
            
            % build frequency lookup table
            obj.buildFreqTable;
            
        end
        
        function satIds = getSatIds(obj, PRNs, constIds)
            % Get indices of the satellites for a given vector of PRNs and
            % constIds.
            % 
            %   satIds = obj.getSatIds(PRNs, constIds)
            %   
            
            % check input dimensions
            if length(PRNs) ~= length(constIds)
                error('Need equal length of PRN and constId vectors.')
            end
            % I acutlaly want to work with row vectors 
            if iscolumn(PRNs)
                PRNs = PRNs';
            end
            if iscolumn(constIds)
                constIds = constIds';
            end
            
            % get sat index mapping
            satIds = find(obj.satPRN == PRNs & obj.satConstId == constIds) ...
                   - obj.numSats * (0:1:length(PRNs)-1)';
            
        end
        
        function [obsData, satIds] = readRinexData(obj, rinexStruct, ep)
            % Parses measurement data from Rinex-like format.
            %   Returns data in struct ready to be processed by this
            %   navigation engine. Accepts inputs for a single or for
            %   multiple epochs. If data for multiple epochs is contained
            %   in the rinexStruct, a second input indicating which epoch
            %   to parse is required (defaults to 1).
            %   
            %   [obsData, satIds] = obj.readRinexData(rinexObsStruct, ep)
            %   
            %   Inputs:
            %   rinexStruct
            %       .meas       - struct of measurements sorted by their
            %                   RINEX 3 codes (C1C, ...) each of size N x M
            %       .PRN        - N x 1 vector of satellite PRNs
            %       .constInds  - N x 1 vector of constellation indices
            %       .tLock      - (optional) N x M carrier phase lock time
            %   Equal input form as to navsu.ppp.preprocessGnssObs()
            %   ep              - (optional) index of the measurement epoch
            %   
            %   Outputs:
            %   obsData         - struct of measurements stored in N x 2
            %                     matricies for N satellites, 2 frequencies
            %       .code       - code phase measurements (m)
            %       .freq       - frequency of each measurement (Hz)
            %       .doppler    - doppler measurement
            %       .carrier    - carrier phase measurement
            %       .CN0        - carrier to noise density ratio (dB-Hz)
            %       .tLock      - carrier phase lock time
            %   satIds          - indices of the N satellites among the
            %                     object's list of satellites
            
            if nargin < 3
                % assume there's only one epoch
                ep = 1;
            end
            
            fn = fieldnames(rinexStruct.meas);
            
            % initialize obs struct
            codeMeas = cellfun(@(x) strcmp(x(1), 'C') && length(x) == 3, fn);
            nSig = length(unique(cellfun(@(x) str2double(x(2)), fn(codeMeas))));
            obsData = struct('code', NaN(length(rinexStruct.PRN), nSig), ...
                             'freq', NaN(length(rinexStruct.PRN), nSig), ...
                             'carrier', NaN(length(rinexStruct.PRN), nSig), ...
                             'tLock', NaN(length(rinexStruct.PRN), nSig), ...
                             'doppler', NaN(length(rinexStruct.PRN), nSig), ...
                             'CN0', NaN(length(rinexStruct.PRN), nSig));
                             
            
            % scan all code meas to analyze each signal
            for fni = find(codeMeas)'
                % get signal identifyer
                sigId = fn{fni}(2:end);
                
                rnxCode = rinexStruct.meas.(fn{fni});
                % for which satellites do I have measurements?
                measIds = find(isfinite(rnxCode(:, ep)) ...
                               & rnxCode(:, ep)~=0);
                
                % how many signals do I have from these satellites already?
                mI = 1;
                while any(isfinite(obsData.code(measIds, mI)))
                    mI = mI + 1;
                end
                
                % now assign the measurements
                obsData.code(measIds, mI) = rnxCode(measIds, ep);
                % assign frequencies
                for mId = 1:length(measIds)
                    fKey = [obj.constLetter(rinexStruct.constInds(measIds(mId))) fn{fni}(2)];
                    obsData.freq(measIds(mId), mI) = obj.freqMap(fKey);
                end
                % assign carrier, doppler, CN0, lock time
                if isfield(rinexStruct.meas, ['L' sigId]) && ~isempty(rinexStruct.meas.(['L' sigId]))
                    obsData.carrier(measIds, mI) = rinexStruct.meas.(['L' sigId])(measIds, ep);
                end
                if isfield(rinexStruct.meas, ['D' sigId]) && ~isempty(rinexStruct.meas.(['D' sigId]))
                    obsData.doppler(measIds, mI) = rinexStruct.meas.(['D' sigId])(measIds, ep);
                end
                if isfield(rinexStruct.meas, ['S' sigId]) && ~isempty(rinexStruct.meas.(['S' sigId]))
                    obsData.CN0(measIds, mI) = rinexStruct.meas.(['S' sigId])(measIds, ep);
                end
                if isfield(rinexStruct, 'tLock') && ~isempty(rinexStruct.tLock)
                    obsData.tLock(measIds, mI) = rinexStruct.tLock(measIds, ep);
                end
            end
            
            % limit to two best frequencies
            [~, Ifreq] = maxk(sum(isfinite(obsData.code), 1), nSig);
            
            obsData = structfun(@(x) x(:, Ifreq([1 2])), obsData, ...
                                'UniformOutput', false);
            
            if nargout > 1
                % now get satIds
                satIds = obj.getSatIds(rinexStruct.PRN, ...
                                       rinexStruct.constInds);
            end
            
        end
        
        function propagateOrbits(obj, satIds, measEpoch)
            % Run orbit propagation for a given measurement epoch.
            %   Corrects measurement epoch for current clock biases and
            %   signal propagation time using current position estimate.
            %   
            %   obj.propagateOrbits(satIds, measEpoch);
            %   
            %   Updates internal memory of satellite data.
            
            % get first guess at transmission epoch
            ttx = obj.transmissionTime(satIds, measEpoch);
            
            % make sure we have rought first estimate of sat positions
            s2update = any(~isfinite(obj.satPos(satIds, :)), 2) ...
                     | abs(obj.internal_satEpoch(satIds) - ttx) > 10;

            if any(s2update)
                obj.updateSatData(satIds(s2update), ttx(s2update));
                % udpate time of transmission
                ttx(s2update) = obj.transmissionTime(satIds(s2update), measEpoch);
            end
            
            % do rough position propagation for theorange estimate
            dtEph = ttx - obj.internal_satEpoch(satIds);
            obj.satPos(satIds, :) = obj.satPos(satIds, :) ...
                               + dtEph .* obj.satVel(satIds, :);
            % recompute time of transmission with now updated theoranges
            ttx = obj.transmissionTime(satIds, measEpoch);
                        
            obj.updateSatData(satIds, ttx);
            
            % account for Sagnac effect (earth rotation during signal
            % flight time)
            obj.compensateSagnacEffect(satIds);
            
        end
        
        function [p, v, t] = pvtSolution(obj, varargin)
            % Perform pvt solution for a given set of satellites for a
            % given epoch for a given set of measurements.
            
            [p, t] = obj.positionSolution(varargin{:});
            
            v = obj.velocitySolution(varargin{:});
            
            
        end
        
        function [pos, tBias, R, varargout] = positionSolution(obj, varargin)
            %Update the position estimate based on new measurements.
            %   Takes as input a new set of measurements. Computes a new
            %   position and time bias solution. Can receive single or
            %   double frequency measurements. If carrier phase
            %   measurements are passed performs carrier smoothing.
            %   
            %   [pos, tBias] = obj.positionSolution( ...
            %       satIds, epoch, freq, obsData)
            %   [pos, tBias, R, prr, P, DOP] = obj.positionSolution( ... )
            
            %   INPUTS (for N satellites, on F frequencies):
            %   satIds      N x 1 vector of integers indicating index of
            %               each measurement
            %   epoch       measurement epoch since first GPS epoch (sec)
            %   obsData     Struct with observables. Has the fields:
            %       freq    N x F vector of signal frequencies (Hz)
            %       code    N x F vector of code phase measurements (meter)
            %       carrier N x F vector of carrier phase measurements
            %               (cycles) OPTIONAL (for carrier smoothing)
            %       tLock   N x F vector of time carrier lock has persisted
            %               (sec) OPTIONAL (for carrier smoothing)
            %   
            %   OUTPUTS (for C involved constellations):
            %   pos         3 x 1 ECEF position solution (m)
            %   tBias       Clock bias in (m) for every involved const.
            %   R           3+C x 3+C Covariance matrix of pos, clock sol.
            %   prr         N x 1 A posteriori pseudorange residuals (m)
            %   P           N x N Residuals information matrix (for chi^2 stat)
            %   DOP         3+C x 3+C Dilution of precision matrix
            
            
            % Step 1: preprocess measurements
            [prMeas, prVar, satIds, epoch] = obj.preprocessMeas(varargin{:});
            % all measurement inputs are now N x 2 matricies [f_SF f_IF]
            
            % Step 2: propagate satellite ephemeris (includes Sagnac
            % effect)
            obj.propagateOrbits(satIds, epoch);
            % exclude satellites where the orbit propagation did not work
            noEph = ~isfinite(obj.theoranges(satIds));
            satIds(noEph) = [];
            prMeas(noEph, :) = [];
            prVar(noEph, :) = [];
            
            % get involved constellations
            consts = unique(obj.satConstId(satIds));
            
            % now solve for the position and clock bias
            % apply elevation mask
            goodEl = obj.satEl(satIds) > obj.elevMask ...
                   | isnan(obj.satEl(satIds));
            nStates = 3 + length(consts);
            
            % check if position solution possible, return if not
            if sum(goodEl) < nStates
                % set NaN outputs, quit
                pos = NaN(3, 1);
                tBias = NaN(length(consts), 1);
                R = NaN(nStates, nStates);
                varargout{1} = NaN(length(varargin{1}), 1);
                varargout{2} = NaN(length(varargin{1}));
                varargout{3} = NaN(nStates, nStates);
                return
            end
            
            loopCounter = 0;
            while sum(goodEl) >= nStates && loopCounter <= 20
                
                % Step 3: correct measurement errors
                [errCorr, SigURE] = obj.UREcorrection(satIds, prMeas, prVar);
                
                % calculate the a-priori pseudo range residual prhat
                prhat = prMeas(:, 1) ...
                      - errCorr ...
                      - obj.theoranges(satIds) ...
                      - obj.tBias(obj.satConstId(satIds)) ...
                      + obj.internal_satClkBias(satIds)*navsu.constants.c;

                % Step 4: do least squares update
                [dx, R, varargout{1:nargout-3}] = obj.doLSupdate( ...
                    satIds(goodEl), prhat(goodEl), SigURE(goodEl));

                obj.position = obj.position + dx(1:3);
                obj.tBias(consts) = obj.tBias(consts) + dx(4:end);
                obj.navCov([1:3, 3+consts'], [1:3, 3+consts']) = R;
                
                % next step depends on size of the update
                normDx = norm(dx);
                if normDx < 1e-6
                    % stop update, we have converged
                    break
                elseif normDx > 1e3
                    % large update, redo orbit propagation with new
                    % theoranges
                    obj.propagateOrbits(satIds, epoch);
                end
                % up the counter
                loopCounter = loopCounter + 1;
                
                % re-apply elevation mask
                goodEl = obj.satEl(satIds) > obj.elevMask ...
                       | isnan(obj.satEl(satIds));

            end
            
            if all(isfinite(obj.position))
                obj.tPosSol = epoch;
            end
            
            % compile outputs [pos, tBias, R, prr, P, DOP]
            if nargout > 0
                pos = obj.position;
                if nargout > 1
                    tBias = obj.tBias(consts);
                    if nargout > 3
                        prr = NaN(length(varargin{1}), 1);
                        [~, LocB] = ismember(satIds(goodEl), varargin{1});
                        prr(LocB) = varargout{1};
                        varargout{1} = prr;
                        if nargout > 4
                            P = NaN(length(varargin{1}));
                            P(LocB, LocB) = varargout{2};
                            varargout{2} = P;
                        end
                    end
                        
                end
            end
        end
        
        function [vel, tRate] = velocitySolution(obj, satIds, epoch, obsData)
            %Compute the user velocity and clock bias rates.
            %   Based on the current position solution.
            
            % check timing
            dt = epoch - obj.tPosSol;
            if abs(dt) > 1
                warning(['Computing velocity with outdated position. ', ...
                    'Position solution ', num2str(dt), ' sec old.'])
            end
            
            % retrieve inputs from structure
            % (choose first freq. measurements for now)
            dopplerMeas = obsData.doppler(:, 1);
            freq = obsData.freq(:, 1);
            
            % exclude satellites without measurements
            sNoMeas = all(~isfinite(dopplerMeas) | ~isfinite(freq), 2);
            satIds(sNoMeas) = [];
            dopplerMeas(sNoMeas, :) = [];
            freq(sNoMeas, :) = [];
            
            % Step 1: propagate orbits (if necessary)
            obj.propagateOrbits(satIds, epoch);
            
            % now limit to satellites with orbit info
            sNoEph = isnan(obj.theoranges(satIds)) ...
                   | any(isnan(obj.satVel(satIds, :)), 2);
            satIds(sNoEph) = [];
            dopplerMeas(sNoEph, :) = [];
            freq(sNoEph, :) = [];
            
            
            % Step 2: perform least squares solution
            % get rate residual
            unitvecs = obj.losVectors(satIds, :) ./ obj.theoranges(satIds);
            rateRes = dot(obj.satVel(satIds, :)', -unitvecs')' ...
                    - dopplerMeas .* navsu.constants.c ./ freq;
            % solve LS nav equation
            rateSol = obj.doLSupdate(satIds, rateRes, ...
                                     obj.codeMeasSigma(satIds, freq));
            
            % store in object properties
            obj.velocity = rateSol(1:3);
            obj.tBiasRate(unique(obj.satConstId(satIds))) = rateSol(4:end);
            
            if nargout > 0
                % compile outputs
                vel = obj.velocity;
                if nargout > 1
                    tRate = rateSol(4:end);
                end
            end
            
        end
        
        function G = Gmatrix(obj, satIds)
            %Computes the geometry matrix.
            %   Takes as input the line of sight vectors in coordinate
            %   system of choice and a vector indicating the constellation
            %   of each vector. Generates matrix in ECEF frame.
            %   
            %   G = Gmatrix(losVectors, constellations)
            %   
            %   Inputs:
            %   satIds      vector indicating which satellites are to be
            %               included
            
            if nargin < 2
                satIds = true(size(obj.satPRN));
            end
            
            % normalize los vectors
            unitvecs = obj.losVectors(satIds, :) ./ obj.theoranges(satIds);
            
            % create geometry matrix
            consts = unique(obj.satConstId(satIds));
            timeMat = zeros(length(unitvecs), numel(consts));
            for c = 1:length(consts)
                timeMat(obj.satConstId(satIds) == consts(c), c) = 1;
            end
            G = [-unitvecs, timeMat];
            
        end
        
        function dop = DOP(obj, satIds)
            %Calculate Dilution of Precision (DOP) matrix.
            
            % limit so satellites with position
            if nargin < 2
                satIds = find(isfinite(obj.theoranges));
            else
                satIds(~isfinite(obj.theoranges(satIds))) = [];
            end
            
            G = obj.Gmatrix(satIds);
            dop = inv(G'*G);
        end
        
        function [errCorr, SigURE] = UREcorrection(obj, satIds, prMeas, measVar)
            %Compute User Range Error (URE) corrections and Variances.
            
            % retrieve TGD correction
            TGD = obj.satTGD(satIds);
            TGD(isnan(TGD)) = 0; % mostly for GLONASS
            
            % relativistic clock correction
            % TODO move this to propagation method to correct sat clock
            % directly? Check IS document to see what is recommended there.
            relCorr = 2 * sum(obj.satPos(satIds, :) ...
                           .* obj.satVel(satIds, :), 2) ...
                           / navsu.constants.c;
                  
            
            % Control system accuracy
            SigCS = obj.satAcc(satIds);
            
            % estimate iono delay
            llh = obj.positionLLH';
            if all(isfinite(llh)) && llh(3) > -1e3
                % have somewhat realistic position, estimate iono and tropo
                el = obj.satEl(satIds);
                az = obj.satAz(satIds);
                if ~isempty(obj.Beph)
                    ionoDelay = klobucharModel(obj.Beph.iono, ...
                                           obj.internal_satEpoch(satIds), ...
                                           llh(1)/180*pi, llh(2)/180*pi, ...
                                           az, el) * navsu.constants.c;
                else
                    ionoDelay = NaN(size(az));
                end
            
                % get tropo error
                if abs(llh(3)) < 1e5
                    doy = navsu.time.jd2doy(navsu.time.epochs2jd( ...
                        obj.internal_satEpoch(satIds)));
                    params.tropModel = 'UNB3';
                    [tropo,~,~] = navsu.ppp.models.tropDelay( ...
                        el*180/pi, az*180/pi, ...
                        llh(:,3), llh(:,1), llh(:,2), doy, params, [], [], []);
                else
                    tropo = NaN(sum(satIds > 0), 1);
                end
            else
                [ionoDelay, tropo] = deal(NaN(sum(satIds > 0), 1));
                
                % take conservative elevation value for propagation delay
                el = 15/180*pi * ones(sum(satIds > 0), 1);
            end
            
            % get Dual Frequency iono delay estimate (this includes TGD!)
            ionoDelayDF = prMeas(:, 1) - prMeas(:, end) - TGD;
            
            % retrieve receiver measurement noise from properties
            % receiver noise & multipath
            SigRNM = measVar(:, 1);
            
            % get measurement weighting based on current position estimate
            % propagation error
            SigTropo = (2.5 ./ sqrt(1 - (cos(el)/1.001).^2)).^2;
            SigTropo(isfinite(tropo)) = SigTropo(isfinite(tropo)) * (0.1)^2;
            % Iono error: default is full iono error
            R_E = navsu.constants.rEarth;
            SigIono = 10^2 * (1 - (R_E * sin(pi/2 - el) ./ (R_E+350000)).^2).^(-1);
            % leverage result of Klobuchar model
            SigIono(isfinite(ionoDelay)) = ionoDelay(isfinite(ionoDelay)).^2;
            % determine when to use DF iono delay estimate
            useDFionoDelay = isfinite(ionoDelayDF) ...
                           & measVar(:, end) < measVar(:, 1) + SigIono;
            % apply DF iono delay estimation where appropriate
            SigIono(useDFionoDelay)     = measVar(useDFionoDelay, end);
            ionoDelay(useDFionoDelay)   = ionoDelayDF(useDFionoDelay);
            SigRNM(useDFionoDelay)      = 0; % accounted for in SigIono

            % sum up errors
            errorMatrix = [relCorr, tropo, ionoDelay, TGD];
            errCorr = sum(errorMatrix, 2, 'omitnan');
%             errCorr(all(isnan(errorMatrix), 2)) = NaN;
            
            if nargout > 1
                SigURE = SigRNM ...
                       + SigTropo ...
                       + SigIono ...
                       + SigCS;
            end
        end
        
%         function saveOutputs(obj)
%             % create some way of saving the computed pvt data
%             
%         end
        
    end
    
    methods (Hidden = true, Access = protected)
        % these are not supposed to be accessed by the user
        
        function sigma = codeMeasSigma(obj, satIds, freq)
            % get expected standard deviation of code measurement as
            % function of elevation. This value matches a low-cost receiver.
            sigma = 3.5 ./ (0.1 + sin(obj.satEl(satIds))) ./ freq * navsu.constants.c;
            % TODO this really should use chip length, not frequency
            sigma(isnan(sigma)) = 2; % conservative default of 2 m
        end
        
        function sigma = carrierMeasSigma(~, freq)
            % get expected standard deviation of carrier phase measurement.
            sigma = 0.05 ./ freq * navsu.constants.c;
        end
        
        function initializeEngine(obj, consts)
            % Initialize several nav engine properties based on
            % constellation structure.
            % 
            % obj.initializeEngine(consts)
            % 
            % Inputs:
            % consts    logical vector indicating which consts to include
            
            % create const struct
            c = num2cell(consts);
            consts = navsu.readfiles.initConstellation(c{1:5});
            
            % initialize several properties
            obj.satPRN              = consts.PRN';
            obj.satConstId          = consts.constInds';
            activeConsts            = unique(consts.constInds);
            obj.numConsts           = numel(activeConsts);
            obj.numSats             = consts.nEnabledSat;
            obj.constellations      = consts.constellations(activeConsts);
            obj.tBias               = zeros(obj.numConsts, 1);
            obj.tBiasRate           = zeros(obj.numConsts, 1);
            obj.navCov              = NaN(3+obj.numConsts);
            obj.rateCov             = NaN(3+obj.numConsts);
            obj.satPos              = NaN(consts.nEnabledSat, 3);
            obj.satVel              = NaN(consts.nEnabledSat, 3);
            obj.satAcc              = NaN(consts.nEnabledSat, 1);
            obj.satTGD              = NaN(consts.nEnabledSat, 1);
            obj.internal_satClkBias = NaN(consts.nEnabledSat, 1);
            obj.internal_satClkRate = NaN(consts.nEnabledSat, 1);
            obj.internal_satEpoch   = NaN(consts.nEnabledSat, 1);
            for f_i = 1:length(obj.CS)
                if f_i == 3
                    tMax = 1800; % code-carrier div. no issue for DF
                else
                    tMax = 100; % to avoid code-carrier divergence
                end
                obj.CS(f_i) = ...
                    navsu.lsNav.CarrierSmoother(consts.nEnabledSat, tMax);
            end
            
        end
        
        function buildFreqTable(obj)
            % Build lookup table of frequency for each signal.
            % Can be accessed by two element string indicating
            % constellation and RINEX3 signal number. The constellation is
            % identified as one of 'GRECJ'.
            
            attachLetter = @(a, l) arrayfun(@(x) [l, num2str(x)], a, ...
                                            'UniformOutput', false);
            
            obj.constLetter = 'GRECJ';
            freqKeys = [attachLetter([1 2 5],       obj.constLetter(1)), ...
                        attachLetter([1 4 2 6 3],   obj.constLetter(2)), ...
                        attachLetter([1 5 7 8 6],   obj.constLetter(3)), ...
                        attachLetter([2 1 5 7 8 6], obj.constLetter(4)), ...
                        attachLetter([1 2 5 6],     obj.constLetter(5))];
            
            freqVals = [1575.42 1227.6 1176.45 ...
                        1602 1600.995 1246 1248.06 1202.025 ...
                        1575.42 1176.45 1207.14 1191.795 1278.75 ...
                        1561.098 1575.42 1176.45 1207.14 1191.795 1268.52 ...
                        1575.42 1227.6 1176.45 1278.75]*1e6;
                
            obj.freqMap = containers.Map(freqKeys, freqVals);
            
        end
        
        
        function [prMeas, prVar, satIds, epoch] = preprocessMeas( ...
                obj, satIds, epoch, obsData)
            %Preprocess the passed inputs and measurements
            %   
            %   Brings inputs to desired sizes and sets defaults for
            %   optional carrier phase inputs and second frequency.
            %   
            %   Then consolidates measurements to a single frequency and
            %   one dual frequency measurement. Performs carrier smoothing
            %   if possible.
            
            % retrieve measurements from input struct and bring to right 
            % dimensions
            if size(obsData.code, 2) == 1
                % single frequency measurement
                codeMeas = [obsData.code, NaN(size(obsData.code))];
            else
                codeMeas = obsData.code;
            end
            
            if size(obsData.freq, 2) == 1
                freq = [obsData.freq, NaN(size(obsData.freq))];
            else
                freq = obsData.freq;
            end
            
            if ~isfield(obsData, 'carrier') || isempty(obsData.carrier)
                % don't have carrier measurement
                carrierMeas = NaN(size(codeMeas));
            elseif size(obsData.carrier, 2) == 1
                carrierMeas = [obsData.carrier, NaN(size(obsData.carrier))];
            else
                carrierMeas = obsData.carrier;
            end
            
            if ~isfield(obsData, 'tLock') || isempty(obsData.tLock)
                % don't have carrier measurement
                lockTime = NaN(size(carrierMeas));
            elseif size(obsData.tLock, 2) == 1
                lockTime = [obsData.tLock, NaN(size(obsData.tLock))];
            else
                lockTime = obsData.tLock;
            end
            
            % Step 0: clean up inputs
            % ensure satellite index is integer index
            if islogical(satIds)
                satIds = find(satIds);
            end
            % limit to satellites with measurement
            noPr = all(~isfinite(codeMeas), 2);
            satIds(noPr) = [];
            freq(noPr, :) = [];
            codeMeas(noPr, :) = [];
            carrierMeas(noPr, :) = [];
            lockTime(noPr, :) = [];
            % convert carrier phase to meter
            carrierMeas = carrierMeas .* navsu.constants.c ./ freq;
            
            % Step 1a: build iono-free combinations
            codeMeas = [codeMeas, obj.buildIFobs(codeMeas, freq)];
            carrierMeas = [carrierMeas, obj.buildIFobs(carrierMeas, freq)];
            lockTime(isnan(lockTime)) = 0;
            lockTime = [lockTime, min(lockTime, [], 2)];
            
%             % could consider doppler as check on cycle slips
%             lCp = [obj.CS(1).lastCarrierPhase(satIds) obj.CS(2).lastCarrierPhase(satIds)];
%             lCpE = [obj.CS(1).tLastSmoothing(satIds) obj.CS(2).tLastSmoothing(satIds)];
%             % get rate of change of carrier phase in cycles/sec
%             carrierDelta = (carrierMeas(:, 1:2) - lCp) ./ (epoch - lCpE) .* freq / navsu.constants.c;
%             carrierVSdoppler = carrierDelta + obsData.doppler(~noPr,  :);
%             % should be close to 0
            
            coSig = obj.codeMeasSigma(satIds, freq).^2;
            coSig = [coSig, obj.buildIFvar(coSig, freq)];
            caSig = obj.carrierMeasSigma(freq).^2;
            caSig = [caSig, obj.buildIFvar(caSig, freq)];
            
            % Step 1b: perform carrier smoothing
            [prSmoothed, prVarRcvr] = deal(NaN(size(codeMeas)));
            for f_i = 1:size(carrierMeas, 2)
                [prSmoothed(:, f_i), prVarRcvr(:, f_i)] = obj.CS(f_i).smoothen( ...
                    satIds, epoch, ...
                    codeMeas(:, f_i), carrierMeas(:, f_i), ...
                    lockTime(:, f_i), coSig(:, f_i), caSig(:, f_i));
            end
            % all measurement inputs are now N x 3 matricies [f_1 f_2 f_IF]
            
            % Step 1c: limit to 1 single frequency (SF) measurement
            chooseSecondSF = isfinite(prSmoothed(:, 2)) ...
                           & (prVarRcvr(:, 2) < prVarRcvr(:, 1) ...
                              | isnan(prSmoothed(:, 1)));
            prSF = prSmoothed(:, 1);
            prSF(chooseSecondSF) = prSmoothed(chooseSecondSF, 2);
            prVarSF = prVarRcvr(:, 1);
            prVarSF(chooseSecondSF) = prVarRcvr(chooseSecondSF, 2);
            
            prMeas = [prSF, prSmoothed(:, end)];
            prVar = [prVarSF, prVarRcvr(:, end)];
            
        end
        
        function obsIF = buildIFobs(~, obs, freq)
            % Calculate the iono-free (IF) combination of a measurement for
            % two given frequencies.
            % 
            % obsIF = obj.buildIFobs(obs, freq)
            % 
            % obs   N x 2 matrix of observations (code or carrier) in meter
            % freq  N x 2 matrix of frequencies in Hz
            
            if size(obs, 2) == 2 && size(freq, 2) == 2
                gamma = (freq(:, 1) ./ freq(:, 2)).^2;
                obsIF = (gamma.*obs(:, 1) - obs(:, 2)) ./ (gamma-1);
            else
                obsIF = NaN(size(obs, 1), 1);
            end
        end
        
        function varIF = buildIFvar(~, measVar, freq)
            % Calculate the variance of the iono-free (IF) combination of 
            % a measurement for two given variances and frequencies.
            % 
            % varIF = obj.buildIFobs(measVar, freq)
            % 
            % measVar   N x 2 matrix of meas variances (code or carrier) in
            %           meter^2
            % freq      N x 2 matrix of frequencies in Hz
            if size(measVar, 2) == 2 && size(freq, 2) == 2
                gamma = (freq(:, 1) ./ freq(:, 2)).^2;
                varIF = (gamma.^2.*measVar(:, 1) + measVar(:, 2)) ./ (gamma-1).^2;
            else
                varIF = NaN(size(measVar, 1), 1);
            end
        end
        
        function [satPosRot, satVelRot] = compensateSagnacEffect(obj, satIds)
            % Compute the satellite positions when rotated to account for
            % the Sagnac effect (earth rotation during signal flight time).
            % Reference: see IS GPS 200, 20.3.3.4.3.3.2
            
            We=7.2921151467e-5; %WGS-84 Earth Rotation Rate (rad/sec)
            
            if nargin < 2
                satIds = true(size(obj.satPRN));
            end
            
            %Rotation angle (radians):
            theta = We * obj.theoranges(satIds) / navsu.constants.c;
            sT = sin(theta);
            cT = cos(theta);
                        
            % do the rotation vectorized for all sats at once
            satPosRot = ...
                [dot([cT, sT, zeros(size(sT))], obj.satPos(satIds, :), 2), ...
                 dot([-sT, cT, zeros(size(sT))], obj.satPos(satIds, :), 2), ...
                 obj.satPos(satIds, 3)];
            
            % also rotate velocity (?)
            satVelRot = ...
                [dot([cT, sT, zeros(size(sT))], obj.satVel(satIds, :), 2), ...
                 dot([-sT, cT, zeros(size(sT))], obj.satVel(satIds, :), 2), ...
                 obj.satVel(satIds, 3)];
             
            % update object property
            obj.satPos(satIds, :) = satPosRot;
            obj.satVel(satIds, :) = satVelRot;
        end
        
        function updateSatData(obj, satIds, epoch)
            %Propagate the selected satellites to chosen epoch.
            % Update internal memory accordingly. Uses precise or broadcast
            % orbit info depending on availability.
            
            % limit to satellites that actually need to be propagated
            n2p = abs(epoch - obj.internal_satEpoch(satIds)) > 1e-9 ...
                | isnan(obj.internal_satEpoch(satIds));
            if ~any(n2p)
                % no need to propagate!
                return
            end
            
            if isempty(obj.Peph)
                % Don't have precise ephemeris. Propagate the broadcast.
                obj.updateSatDataFromBrdc(satIds(n2p), epoch(n2p));
                
            else
                % Propagate the precise ephemeris
                % need: satPos, satVel, satClk.bias, satClk.drift, posVar,
                % TGD
                obj.updateSatDataFromPProd(satIds(n2p), epoch(n2p));
            end
            
        end
        
        function updateSatDataFromBrdc(obj, s2update, epoch)
            %Propagate the navigation broadcast for selected satellites to
            %chosen epoch. Update internal memory accordingly.
            
            % here I could check the memory to see what I need to propagate                
                
            % propagate navigation broadcast to desired epoch
            satInfo = navsu.geo.propNavMsg(obj.Beph, ...
                                           obj.satPRN(s2update), ...
                                           obj.satConstId(s2update), ...
                                           epoch);
            % save everything to internal memory
            obj.satPos(s2update, :) = [satInfo.x, ...
                                       satInfo.y, ...
                                       satInfo.z];
            obj.satVel(s2update, :) = [satInfo.x_dot, ...
                                       satInfo.y_dot, ...
                                       satInfo.z_dot];
            obj.internal_satEpoch(s2update)   = epoch;
            obj.internal_satClkBias(s2update) = satInfo.clock_bias;
            obj.internal_satClkRate(s2update) = satInfo.clock_drift;
            % CS is better than it promises:
            obj.satAcc(s2update)              = satInfo.accuracy.^2 / 4;
            obj.satAcc(s2update(obj.satConstId(s2update) == 2)) = 9; % GLO
            obj.satTGD(s2update) = satInfo.TGD * navsu.constants.c; % in meter
            
        end
        
        function updateSatDataFromPProd(obj, s2update, epoch)
            %Propagate the precise orbit and clock for selected satellites 
            %to chosen epoch. Update internal memory accordingly.
            %   
            %   updateSatDataFromPProd(obj, s2update, epoch)
            %   
            %   Inputs:
            %   s2update
            
            N = sum(s2update>0);
            [obj.satPos(s2update, :), obj.satVel(s2update, :)] = ...
                obj.Peph.propagate(obj.satPRN(s2update), ...
                                   obj.satConstId(s2update), ...
                                   epoch.*ones(N, 1));
            svClock = obj.Peph.clock(obj.satPRN(s2update), ...
                                     obj.satConstId(s2update), ...
                                     epoch.*ones(N, 1));

            obj.internal_satClkRate(s2update) = ...
                (svClock - obj.internal_satClkBias(s2update)) ...
                ./ (epoch - obj.internal_satEpoch(s2update));

            obj.internal_satEpoch(s2update)   = epoch;
            obj.internal_satClkBias(s2update) = svClock;
            % set some high accuracy for orbits
            obj.satAcc(s2update)              = 0.1;

            obj.satTGD(s2update) = 0; % could use DCBs?!
            
        end
        
        function ttx = transmissionTime(obj, satIds, measEpoch)
            %Compute time of signal transmission.
            constIds = obj.satConstId(satIds);
            rcvrClkBias = obj.tBias(constIds) ...
                + (measEpoch - obj.tPosSol) .* obj.tBiasRate(constIds);
            
            ttx = sum([measEpoch*ones(sum(satIds > 0), 1), ...
                - rcvrClkBias / navsu.constants.c, ...
                - obj.theoranges(satIds) / navsu.constants.c],  2, 'omitnan');
            
        end
        
        function [dx, Rout, prr, P, DOP] = doLSupdate(obj, satIds, prhat, measVar)
            % Solve the weighted least squares navigation equation for a
            % given a posteriori residual and measurement variance.
            % 
            %   Applies an elevation mask angle.
            %   
            %   [dx, R] = obj.doLSupdate(satIds, prhat, measVar)
            %   [dx, R, prr, P, DOP] = obj.doLSupdate(satIds, prhat, measVar)
            
            % number of states of the outputs
            outConsts = unique(obj.satConstId(satIds));
            nOutStates = 3 + length(outConsts);
            
            % apply elevation mask
            measExclude = obj.satEl(satIds) < obj.elevMask ...
                        | isnan(prhat) ...
                        | isnan(measVar);
            satIds(measExclude) = [];
            prhat(measExclude) = [];
            measVar(measExclude) = [];
            
            usedConsts = unique(obj.satConstId(satIds))';
            outStates = [1:3, 3+find(ismember(usedConsts, outConsts))];
            
            % initialize outputs
            N = length(measExclude);
            dx = NaN(nOutStates, 1);
            Rout = NaN(nOutStates);
            prr = NaN(N, 1);
            P = NaN(N);
            DOP = NaN(nOutStates, nOutStates);
            
            % solve weighted least squares equation
            G = obj.Gmatrix(satIds); % geometry matrix
            W = diag(1./measVar); % weighting matrix
            [U,S,V] = svd(G'*W*G);
            
            % check if it can be solved
            if S(1, 1) / S(end, end) > 1e8
                %geometry too poorly conditioned. Don't do an update.
                warning('Deficient geometry, position update impossible.')
                return
            end

            % a posteriori position covariance matrix R = inv(G'*W*G)
            R = V * diag(1./diag(S)) * U';
            LSest = R * G' * W; % LS estimator
            dx(outStates) = LSest * prhat;
            
            % calculate optional outputs
            if nargout > 1
                Rout(outStates, outStates) = R;
            end
            if nargout > 2
                % a-posteriori pseudo range residual
                prr(~measExclude) = prhat - G*dx;
                if nargout > 3
                    % residual information matrix
                    P(~measExclude, ~measExclude) = W - W * G * LSest;
                    if nargout > 4
                        % dilution of precision
                        DOP(outStates, outStates) = inv(G'*G);
                    end
                end
            end
        end
        
    end
end

