classdef DFMCnavigationEngine < matlab.mixin.Copyable
    %DFMCnavigationEngine Dual frequency, multi constellation navigation.
    %   Computes navigation solutions from dual frequency, multi
    %   constellation GNSS measurements. Can compute position, velocity and
    %   time solutions.
    %   Offers options to use/not use models for iono and tropo delay as
    %   well as carrier smoothing. Can be run with precise or broadcast
    %   orbits. Options are being set by obj.use... properties.
    %   
    %   Properties:
    %   velocity            current velocity estimate in ECEF in (m/s)
    %   tBias               clock bias for each constellation in (s)
    %   tBiasRate           clock bias rate for each constellation in (s/s)
    %   tPosSol             time of last position solution in sec since first GPS epoch
    %   usePreciseProducts  should precise orbit products be used?
    %   useCarrierSmoothing attempt carrier smoothing?
    %   useIonoModel        use Klobuchar model to estimate iono delay?
    %   useTropoModel       use UNB3 model to estimate tropo delay
    %   smoothingConstantIF iono-free carrier smoothing time constant in sec
    %   smoothingConstant   single frequency carrier smoothing time constant in sec
    %   numConsts           number of included constellations
    %   activeConsts        indices of the involved constellations
    %   numSats             number of included satellites
    %   constellations      char vector indicating the useable constellations
    %   navCov              Position, clock bias covariance matrix
    %   rateCov             Velocity, clock rate covariance matrix
    %   theoranges          ranges user to satellite (length of losVectors)
    %   satAcc              Variance of satellite position
    %   satPRN              PRN of each satellite
    %   satConstId          constellation id of each satellite
    %   satTGD              Timing group delay of each satellite in sec 
    %                       (see IS ﻿20.3.3.3.3.2)
    %   satEl               elevation of satellites in deg
    %   satAz               azimuth of satellites in deg
    %   CS                  array of objects for carrier smoothing.  One
    %                       smoother for each signal, one for each Dual Freq.
    %   elevMask            elevation mask angle in rad
    %   position            current position estimate in ECEF coordinates in (m)
    %   positionLLH         position estimate in lat, lon, height (deg, deg, m)
    %   losVectors          line of sight vectors user to satellite
    %   satPos              satellite position at time of broadcast
    %   satVel              satellite velocity at time of broadcast
    %   satEph              svOrbitClk object containing ephemeris data
    
    
    properties
        velocity (3,1)   % current velocity estimate in ECEF in (m/s)
        tBias      % time bias for each constellation in (s)
        tBiasRate   % time bias rate for each constellation in (s/s)
        tPosSol = NaN; % time of last position solution in sec since first GPS epoch
        usePreciseProducts (1,1) logical = false;
        useCarrierSmoothing (1,1) logical = true;
        useIonoModel (1,1) logical = true;
        useTropoModel (1,1) logical = true;
        smoothingConstantIF (1,1) double {mustBeReal, mustBeFinite} = 1800;
        smoothingConstant (1,1) double {mustBeReal, mustBeFinite} = 100;
        numConsts   % number of included constellations
        activeConsts% indices of the involved constellations
        numSats     % number of included satellites
        constellations % char vector indicating the useable constellations
        navCov      % Position, time bias covariance matrix
        rateCov     % Velocity, time rate covariance matrix
        theoranges  % ranges user to satellite (length of losVectors)
        satAcc      % Variance of satellite position
        satPRN      % PRN of each satellite
        satConstId  % constellation id of each satellite
        satTGD      % Timing group delay of each satellite in sec (see IS ﻿20.3.3.3.3.2)
        satEl       % elevation of satellites in deg
        satAz       % azimuth of satellites in deg
        CS navsu.lsNav.CarrierSmoother % array of objects for carrier 
        % smoothing. One smoother for each signal, one for each Dual Freq.
        elevMask = 15*pi/180;    % elevation mask angle
    end
    
    properties (Dependent)
        position    % current position estimate in ECEF coordinates in (m)
        positionLLH % position estimate in lat, lon, height (deg, deg, m)
        losVectors  % line of sight vectors user to satellite
        satPos      % satellite position at time of broadcast
        satVel      % satellite velocity at time of broadcast
        satEph      % svOrbitClk object containing ephemeris data
    end
    
    properties (Access = private, Hidden = true)
        % internal memory
        internal_position (3,1)
        internal_losVectors
        internal_satPos = zeros(0, 3);
        internal_satVel = zeros(0, 3);
        internal_satClkBias
        internal_satClkRate
        internal_satEpoch       % epoch of last orbit propagation
        internal_svOrbClk navsu.svOrbitClock % satellite ephemeris data
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
        
        function eph = get.satEph(obj)
            % retrieve from memory
            eph = obj.internal_svOrbClk;
        end
        function set.satEph(obj, ephData)
            % Save the passed satellite ephemeris information. Detect if
            % it's an eph struct or a svOrbitClock object and handle
            % accordingly.

            if isa(ephData, 'navsu.svOrbitClock')

                % got svOrbitClock product directly
                obj.internal_svOrbClk = ephData;

            elseif isstruct(ephData)
                % eph struct was passed
                % first retrieve for which constellations is has data
                useConst = false(1, 5);
                % check each constellation to see if it can be used
                for cI = 1:length(useConst)
                    cStr = lower(navsu.svprn.convertConstIndName(cI));
                    useConst(cI) = isfield(ephData, cStr) ...
                                && ~isempty(ephData.(cStr)) ...
                                && numel(fieldnames(ephData.(cStr))) > 1;
                end
                % initialize and populate svOrbitClock object
                obj.internal_svOrbClk = navsu.svOrbitClock('constUse', ...
                                                           useConst);
                obj.internal_svOrbClk.BEph = ephData;
                
            end

            % also initialize internal sat memory properties accordingly
            obj.initializeEngine;
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
        function obj = DFMCnavigationEngine(ephData, x0)
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
                % store ephemeris info
                obj.satEph = ephData;
            end
            
            % read further inputs
            if nargin > 1
                obj.position = x0;
            else
                obj.position = zeros(3, 1);
            end
            
        end
        
        function [pos, vel, tBias, R, prr, P, chi2stat, dop, useDF] = ...
                batchPvtSolution(obj, obsGnss, useConst, fBands)
            % Run the PVT solver for measurements at multiple epochs.
            %   Automatically runs the recursive least squares positioning
            %   over a batch of measurements at multiple epochs.
            %   
            %   [posECEF, velECEF, tBias, R, prr, P, chi2stat, dop] = ...
            %       obj.batchPvtSolution(obsGnss)
            %   
            %   Can be called with a limited number of constellations or
            %   frequencies:
            %   [pos, vel, tBias, R, prr, P, chi2s, dop] = ...
            %       obj.batchPvtSolution(obsGnss, useConst, freqs)
            %   
            %   Inputs:
            %   obsGnss     Struct of measurement data, same input as for
            %               obj.readRinexData and 
            %               navsu.ppp.preprocessGnssObs. Contains data for
            %               N satellites at M epochs.
            %       .constInds  N x 1 vector of indices indicating the
            %                   constellation of each signal
            %       .PRN        N x 1 vector of satellite PRNs
            %       .epochs     1 x M vector of measurement epochs
            %       .meas       struct of measurements sorted by their
            %                   RINEX 3 codes (C1C, ...) each of size N x M
            %       .tLock      (optional) struct of carrier phase lock
            %                   time, same fields as .meas struct.
            %   useConst    OPTIONAL indices limiting the constellations to
            %               be used. Default is 1:5.
            %   fBands      OPTIONAL indices indicating the frequency bands
            %               to be used. Default is [1 2].
            %   
            %   Outputs:
            %   pos     3 x M matrix of ECEF position solutions
            %   vel     3 x M matrix of ECEF velocity solutions
            %   tBias   C x M matrix of clock biases for C constellations
            %   R       (3+C) x (3+C) x M position, time covariance matrix
            %   prr     N x M matrix of pseudorange residuals
            %   P       N x N x M matrix of residual information matricies
            %   chi2s   1 x N vector of residual chi^2 statistic values
            %   dop     (3+C) x M matrix of DOP values (diagonal values of
            %           DOP matricies)
            %   useDF   1 x M logical vector that is positive where a dual
            %           frequency position solution was computed.

            if nargin < 3 || isempty(useConst)
                % default to using everything
                useConst = 1:5;
            elseif islogical(useConst) || all(ismember(useConst, [0 1]))
                % ensure integer indices
                useConst = find(useConst);
            end
            % limit to constellations in the data
            useConst = intersect(useConst, obsGnss.constInds);

            if nargin < 4 || isempty(fBands)
                % default to using dual frequency
                fBands = 1:2;
            elseif islogical(fBands)
                % ensure integer indices
                fBands = find(fBands);
            end

            % select the satellites to use
            sats = ismember(obsGnss.constInds, useConst) ...
                 & any(obsGnss.constInds == obj.satConstId' ...
                       & obsGnss.PRN == obj.satPRN', 2);
            
            % initialize all the outputs
            nEpochs = length(obsGnss.epochs);
            nConsts = length(useConst);
            nSats = sum(sats);

            pos = NaN(3, nEpochs);
            vel = NaN(3, nEpochs);
            tBias = NaN(nConsts, nEpochs);
            prr = NaN(nSats, nEpochs); % pseudorange residuals
            R = NaN(3+nConsts, 3+nConsts, nEpochs); % cov. matrix
            P = NaN(nSats, nSats, nEpochs); % residuals information matrix
            chi2stat = NaN(1, nEpochs); % chi-squared statistic
            dop = NaN(3+nConsts, nEpochs); % geometric dilution of precision
            useDF = false(1, nEpochs); % are we using the DF solution?

            for ep = 1:nEpochs
                % parse obs data for this epoch, retrieve satellite indices
                [obsData, satIds] = obj.readRinexData(obsGnss, sats, ep);

                % limit to selected frequencies
                obsData = obj.frequencyMask(obsData, fBands);
                
                % do PVT computation
                [pos(:, ep), tBias(:, ep), R(:, :, ep), prr(:, ep), ...
                    P(:, :, ep), dopMat, useDF(ep)] = ...
                    obj.positionSolution(satIds, obsGnss.epochs(ep), obsData);
                
                vel(:, ep) = ...
                    obj.velocitySolution(satIds, obsGnss.epochs(ep), obsData);
                
                % also compute chi square statistic, save dop
                s = isfinite(prr(:, ep));
                chi2stat(ep) = prr(s, ep)' * P(s, s, ep) * prr(s, ep);
                dop(:, ep) = diag(dopMat);
            end

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
            PRNmatches = obj.satPRN == PRNs & obj.satConstId == constIds;
            satIds = find(PRNmatches) ...
                   - obj.numSats * (find(any(PRNmatches, 1))-1)';
            
        end
        
        function [obsData, satIds] = readRinexData(obj, rnxStruct, sats, ep)
            % Parses measurement data from Rinex-like format.
            %   Returns data in struct ready to be processed by this
            %   navigation engine. Accepts inputs for a single or for
            %   multiple epochs. If data for multiple epochs is contained
            %   in the rinexStruct, a second input indicating which epoch
            %   to parse is required (defaults to 1).
            %   
            %   [obsData, satIds] = obj.readRinexData(rinexObsStruct, sats, ep)
            %   
            %   Inputs:
            %   rinexStruct
            %       .meas       - struct of measurements sorted by their
            %                   RINEX 3 codes (C1C, ...) each of size N x M
            %       .PRN        - N x 1 vector of satellite PRNs
            %       .constInds  - N x 1 vector of constellation indices
            %       .tLock      - (optional) struct of carrier phase lock
            %                   time, fields as in .meas struct that start
            %                   with 'L'.
            %   Equal input form as to navsu.ppp.preprocessGnssObs()
            %   sats            - (optional) indices of satellites to use
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
            
            if nargin < 3 || isempty(sats)
                sats = 1:length(rnxStruct.PRN);
            elseif islogical(sats)
                % ensure integer indices
                sats = find(sats);
            end
            if nargin < 4
                % assume there's only one epoch
                ep = 1;
            end

            % limit to satellites handled by this object
            validSats = any(obj.satPRN == rnxStruct.PRN(sats)' ...
                          & obj.satConstId == rnxStruct.constInds(sats)', 1);
            sats = sats(validSats);
            
            % identify code measurements among observables
            [codeObservables, fn] = obj.findRnxCodeObs(rnxStruct, sats, ep);

            % now we can preallocate the output struct
            nSat = length(sats);
            nSig = max(sum(codeObservables, 1));
            obsData = struct('code',    NaN(nSat, nSig), ...
                             'freq',    NaN(nSat, nSig), ...
                             'carrier', NaN(nSat, nSig), ...
                             'tLock',   NaN(nSat, nSig), ...
                             'doppler', NaN(nSat, nSig), ...
                             'CN0',     NaN(nSat, nSig), ...
                             'fBand',   NaN(nSat, nSig));
            obsData.rnxCode = repmat({''}, nSat, nSig);
            
            for c = 1:size(codeObservables, 2)
                % get signal names of code measurements for this constellation
                % do so in increasing number in case some have to be cut
                % off
                fni = sort(fn(codeObservables(:, c)));
            
                % scan all code meas to analyze each signal
                for sI = 1:length(fni)
                    % get signal identifyer
                    codeId = fni{sI};
                    sigId = codeId(2:end);
                    
                    rnxCodeMeas = rnxStruct.meas.(codeId)(sats, ep);
                    % for which satellites do I have measurements?
                    satIds = find(isfinite(rnxCodeMeas) & rnxCodeMeas ~= 0);

                    if isempty(satIds)
                        % no measurements on this frequency
                        continue;
                    end
                    
                    % now assign the measurements
                    obsData.code(satIds, sI) = rnxCodeMeas(satIds);
                    sIds = sats(satIds);

                    % assign frequencies
                    obsData.freq(satIds, sI) = ...
                        navsu.svprn.mapSignalFreq(str2double(codeId(2)), ...
                                                  rnxStruct.PRN(sIds), ...
                                                  rnxStruct.constInds(sIds), ...
                                                  navsu.time.epochs2jd(rnxStruct.epochs(ep)));
                    % assign carrier, doppler, CN0, carrier lock time
                    if isfield(rnxStruct.meas, ['L' sigId]) ...
                            && ~isempty(rnxStruct.meas.(['L' sigId]))
                        obsData.carrier(satIds, sI) = rnxStruct.meas.(['L' sigId])(sIds, ep);
                    end
                    if isfield(rnxStruct.meas, ['D' sigId]) ...
                            && ~isempty(rnxStruct.meas.(['D' sigId]))
                        obsData.doppler(satIds, sI) = rnxStruct.meas.(['D' sigId])(sIds, ep);
                    end
                    if isfield(rnxStruct.meas, ['S' sigId]) ...
                            && ~isempty(rnxStruct.meas.(['S' sigId]))
                        obsData.CN0(satIds, sI) = rnxStruct.meas.(['S' sigId])(sIds, ep);
                    end
                    if isfield(rnxStruct, 'tLock') ...
                            && ~isempty(rnxStruct.tLock) ...
                            && isfield(rnxStruct.tLock, ['L' sigId]) ...
                            && ~isempty(rnxStruct.tLock.(['L' sigId]))
                        obsData.tLock(satIds, sI) = rnxStruct.tLock.(['L' sigId])(sIds, ep);
                    end
                    % store signal identifier
                    obsData.rnxCode(satIds, sI) = {sigId};
                    % store frequency band
                    obsData.fBand(satIds, sI) = str2double(sigId(1));
                end
            end
            
            if nargout > 1
                % now get satIds
                satIds = obj.getSatIds(rnxStruct.PRN(sats), ...
                                       rnxStruct.constInds(sats));
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
        
        function [pos, tBias, R, varargout] = positionSolution(obj, ...
                allSatIds, epoch, obsData)
            %Update the position estimate based on new measurements.
            %   Takes as input a new set of measurements. Computes a new
            %   position and time bias solution. Can receive single or
            %   double frequency measurements. If carrier phase
            %   measurements are passed performs carrier smoothing.
            %   
            %   [pos, tBias] = obj.positionSolution( ...
            %       satIds, epoch, obsData)
            %   [pos, tBias, R, prr, P, DOP, useDF] = obj.positionSolution( ... )
            
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
            %   useDF       logical indicating the use of dual frequency
            
            
            % Step 1: preprocess measurements
            [prMeas, prVar, satIds, freqs] = obj.preprocessMeas( ...
                allSatIds, epoch, obsData);
            % all measurement inputs are now N x 2 matricies [f_SF f_IF]

            % get involved constellations
            consts0 = unique(obj.satConstId(allSatIds));

            % update the current estimate of the clock bias
            obj.updateClkBias(epoch);
            
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
            satConstIds = obj.getActiveConstId(obj.satConstId(satIds));
            
            % now solve for the position and clock bias
            % apply elevation mask
            goodEl = obj.satEl(satIds) > obj.elevMask ...
                   | isnan(obj.satEl(satIds));
            nStates = 3 + length(consts);
            
            % check if position solution possible, return if not
            if sum(goodEl) < nStates
                % set NaN outputs, quit
                nStates0 = 3 + length(consts0);
                pos = NaN(3, 1);
                tBias = NaN(length(consts0), 1);
                R = NaN(nStates0, nStates0);
                nSats = length(allSatIds);
                varargout{1} = NaN(nSats, 1);  % prr
                varargout{2} = NaN(nSats); % P
                varargout{3} = NaN(nStates0, nStates0); % dop
                varargout{4} = false; % useDF
                return
            end

            % Step 3: correct measurement errors
            [errCorr, SigURE, varargout{4}] = obj.UREcorrection(...
                satIds, prMeas, prVar, freqs);
            
            
            % initialize loop counter, last state update
            loopCounter = 0; dx = inf;

            while sum(goodEl) >= nStates && loopCounter <= 20
                
                % use latest elevation mask
                goodElNow = goodEl;

                % calculate the a-priori pseudo range residual prhat
                prhat = prMeas(goodElNow, 1) ...
                      - errCorr(goodElNow) ...
                      - obj.theoranges(satIds(goodElNow)) ...
                      - obj.tBias(satConstIds(goodElNow)) ...
                      + obj.internal_satClkBias(satIds(goodElNow))*navsu.constants.c;

                % Step 4: do least squares update
                [dxi, Ri, varargout{1:min(nargout-3, 3)}] = obj.doLSupdate( ...
                    satIds(goodElNow), prhat, SigURE(goodElNow));
                
                % check for oscillating state, bad update
                if all(dxi + dx < 1e-6) || any(isnan(dxi(1:3)))
                    break;
                else
                    dx = dxi;
                end

                obj.position = obj.position + dx(1:3);

                % constellation-depending updates:
                biasUpdate = dx(4:end);
                % check for which consts have been updated
                updatedConsts = obj.getActiveConstId(unique( ...
                    obj.satConstId(satIds(goodElNow))));
                % only update the ones where we have an update
                haveBiasUpdate = isfinite(biasUpdate);
                updatedConsts = updatedConsts(haveBiasUpdate);

                obj.tBias(updatedConsts) = obj.tBias(updatedConsts) ...
                                         + biasUpdate(haveBiasUpdate);
                % same for covariance
                stateVars = [1:3, 3+updatedConsts'];
                compStates = isfinite(dx);
                obj.navCov(stateVars, stateVars) = Ri(compStates, compStates);
                
                % next step depends on size of the update
                normDx = sqrt(sum(dx.^2, 'omitnan'));
                if normDx < 1e-4
                    % stop update, we have converged
                    break;
                elseif normDx > 1e3
                    % large update, redo orbit propagation with new
                    % theoranges
                    obj.propagateOrbits(satIds, epoch);
                    % redo pseudorange error estimates
                    [errCorr, SigURE, varargout{4}] = obj.UREcorrection(...
                        satIds, prMeas, prVar, freqs);
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
                    % tBias
                    tBias = obj.tBias(obj.getActiveConstId(consts0));
                    if nargout > 2
                        % R
                        nStates = 3+length(tBias);
                        R = NaN(nStates);
                        R(stateVars, stateVars) = obj.navCov(stateVars, stateVars);

                        if nargout > 3
                            nSats = length(allSatIds);
                            % make sure prr, P are for the right satIds
                            prr = NaN(nSats, 1);
                            [~, LocB] = ismember(satIds(goodElNow), allSatIds);
                            prr(LocB) = varargout{1};
                            varargout{1} = prr;
                            if nargout > 4
                                P = NaN(nSats);
                                P(LocB, LocB) = varargout{2};
                                varargout{2} = P;
                                if nargout > 5
                                    % dop
                                    tempDop = varargout{3};
                                    varargout{3} = NaN(nStates);
                                    varargout{3}(stateVars, stateVars) = tempDop(compStates, compStates);
                                end
                            end
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
            if dt > 0 || any(~isfinite(obj.satVel(satIds, :)), 'all')
                obj.propagateOrbits(satIds, epoch);
            end
            
            % now limit to satellites with orbit info
            sNoEph = isnan(obj.theoranges(satIds)) ...
                   | any(isnan(obj.satVel(satIds, :)), 2);
            satIds(sNoEph) = [];
            dopplerMeas(sNoEph, :) = [];
            freq(sNoEph, :) = [];
            
            satConstIds = obj.getActiveConstId(obj.satConstId(satIds));
            
            % Step 2: perform least squares solution
            % get rate residual
            unitvecs = obj.losVectors(satIds, :) ./ obj.theoranges(satIds);
            rateRes = dot(obj.satVel(satIds, :)', -unitvecs')' ...
                    - dopplerMeas .* navsu.constants.c ./ freq ...
                    - obj.tBiasRate(satConstIds);
            % solve LS nav equation
            rateSol = obj.doLSupdate(satIds, rateRes, ...
                                     obj.codeMeasSigma(satIds, freq));
            
            % store in object properties
            obj.velocity = rateSol(1:3);
            satConstIdsU = unique(satConstIds);
            obj.tBiasRate(satConstIdsU) = obj.tBiasRate(satConstIdsU) ...
                                        + rateSol(4:end);
            
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
        
        function [errCorr, SigURE, useDF] = UREcorrection(obj, satIds, prMeas, measVar, freqs)
            %Compute User Range Error (URE) corrections and Variances.
            % WARNING: This blindly uses the stored group delay values if a
            % single frequency solution is computed. This can lead to
            % erroneous results when e.g. the wrong Galileo BGD is stored
            % w.r.t. the E5 signal being used.
            
            if nargin < 5
                freqs = 1.57542e9 * ones(size(prMeas, 1), 1);
            end

            gamma = (1.57542e9 ./ freqs(:, 1)).^2;

            % retrieve TGD correction
            TGD = gamma .* obj.satTGD(satIds);
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
                satEpochs = obj.internal_satEpoch(satIds);

                % estimate iono delay using Klobuchar model
                if obj.useIonoModel && ~isempty(obj.satEph.BEph)
                    
                    ionoCorrCoeffs = obj.getIonoCoeffs(mean(satEpochs));
                    
                    ionoDelay = gamma .* navsu.ppp.models.klobuchar( ...
                        ionoCorrCoeffs, satEpochs, ...
                        llh(1)/180*pi, llh(2)/180*pi, ...
                        az, el) * navsu.constants.c;
                else
                    ionoDelay = NaN(size(az));
                end
            
                % get tropo error
                if obj.useTropoModel && abs(llh(3)) < 1e5
                    doy = navsu.time.jd2doy(navsu.time.epochs2jd(satEpochs));
                    params.tropModel = 'UNB3';
                    [tropo, ~, ~] = navsu.ppp.models.tropDelay( ...
                        el*180/pi, az*180/pi, ...
                        llh(:,3), llh(:,1), llh(:,2), doy, params);
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

            haveDF = isfinite(ionoDelayDF) & isfinite(measVar(:, end));
            % use either DF or SF. Use DF if it results in better G'WG
            if sum(haveDF) < 3 + length(unique(obj.satConstId(satIds)))
                useDF = false;
            else
                % compare using SF vs. DF, choose option with smaller
                % largest dimension of R = inv(G'WG): biggest min(eig(G'WG))
                WSF = diag(1./(SigRNM + SigIono + SigTropo + SigCS));
                GSF = obj.Gmatrix(satIds);
                WDF = diag(1 ./ (measVar(haveDF, end) ...
                                 + SigTropo(haveDF) ...
                                 + SigCS(haveDF)));
                GDF = obj.Gmatrix(satIds(haveDF));
                useDF = min(eig(GDF'*WDF*GDF)) > min(eig(GSF'*WSF*GSF));
            end

            if useDF
                % use only DF measurements
                useDFionoDelay = haveDF;
                SigIono(~useDFionoDelay) = inf;
            else
                % don't use DF measurements
                useDFionoDelay = false(size(ionoDelayDF));
            end

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
        
        function initializeEngine(obj)
            % Initialize several nav engine properties based on
            % the involved constellations. Based on the ephemeris
            % information provided to the navigation engine.
            % 
            % obj.initializeEngine
            % 
            % Inputs:
            % none. Uses obj.satEph.
            
            % create const struct
            c = num2cell(obj.satEph.settings.constUse);
            consts = navsu.readfiles.initConstellation(c{1:5});
            
            % initialize several properties
            obj.satPRN              = consts.PRN';
            obj.satConstId          = consts.constInds';
            obj.activeConsts        = unique(obj.satConstId);
            obj.numConsts           = numel(obj.activeConsts);
            obj.numSats             = consts.nEnabledSat;
            obj.constellations      = consts.constellations(obj.activeConsts);
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

            % also set a logical indicator whether to use precise products
            obj.usePreciseProducts = ~isempty(obj.satEph.PEph);

            % make sure the orbit product knows it's using precise data
            if obj.usePreciseProducts
                obj.satEph.orbMode = 'PRECISE';
                if isempty(obj.satEph.PClock)
                    % have to use the broadcast clock
                    obj.satEph.clkMode = 'BROADCAST';
                else
                    obj.satEph.clkMode = 'PRECISE';
                end
            end
            
        end
        
        
        function [prMeas, prVar, satIds, freqs] = preprocessMeas( ...
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
                obsData.code = [obsData.code, NaN(size(obsData.code))];
            end
            
            if size(obsData.freq, 2) == 1
                obsData.freq = [obsData.freq, NaN(size(obsData.freq))];
            end
            
            if ~isfield(obsData, 'carrier') || isempty(obsData.carrier)
                % don't have carrier measurement
                obsData.carrier = NaN(size(obsData.code));
            elseif size(obsData.carrier, 2) == 1
                obsData.carrier = [obsData.carrier, NaN(size(obsData.carrier))];
            end
            
            if ~isfield(obsData, 'tLock') || isempty(obsData.tLock)
                % don't have carrier measurement
                obsData.tLock = zeros(size(obsData.carrier));
            elseif size(obsData.tLock, 2) == 1
                obsData.tLock = [obsData.tLock, zeros(size(obsData.tLock))];
            end
            obsData.tLock(isnan(obsData.tLock)) = 0;

            if ~isfield(obsData, 'fBand') || isempty(obsData.fBand)
                % don't have frequency bands assigned
                obsData.fBand = (1:size(obsData.code, 2)) ...
                                .* ones(size(obsData.code));
                obsData.fBand(~isfinite(obsData.code)) = NaN;
            elseif size(obsData.fBand, 2) == 1
                obsData.fBand = [obsData.fBand, zeros(size(obsData.fBand))];
            end

            if ~isfield(obsData, 'rnxCode') || isempty(obsData.rnxCode)
                % make up rinex signal identifiers based on frequency bands
                fVals = isfinite(obsData.fBand);
                obsData.rnxCode = repmat({''}, size(obsData.fBand));
                obsData.rnxCode(fVals) = arrayfun(@(x) [num2str(x), 'C'], ...
                    obsData.fBand(fVals), 'UniformOutput', false);
            elseif size(obsData.rnxCode, 2) == 1
                obsData.rnxCode = [obsData.rnxCode, ...
                                   repmat({''}, size(obsData.rnxCode))];
            end
            
            % Step 0: clean up inputs
            % ensure satellite index is integer index
            if islogical(satIds)
                satIds = find(satIds);
            end
            % limit to satellites with any measurement
            noPr = all(~isfinite(obsData.code), 2);
            satIds(noPr) = [];
            for fn = fieldnames(obsData)'
                obsData.(fn{1})(noPr, :) = [];
            end

            % convert carrier phase to meter
            obsData.carrier = obsData.carrier ...
                              .* navsu.constants.c ./ obsData.freq;

            % Step 1a: build iono-free combinations
            [obsData, SigRNM] = obj.buildIFobs(satIds, obsData);
            
            % Step 1b: perform carrier smoothing
            [prSmoothed, prVarRcvr] = obj.carrierSmoothing(satIds, epoch, obsData, SigRNM);

            % now the signals could be averaged on each frequency band
%             fBs = unique(obsData.fBand(isfinite(obsData.fBand)));
%             for fbi = 1:length(fBs)
% 
%                 haveSignal = obsData.fBand == fBs(fbi);
%                 nSignals = sum(haveSignal, 2);
%                 prAverage(:, fbi) = sum(prSmoothed .* (haveSignal), 2) ...
%                                     ./ nSignals;
% 
%                 prAvVar(:, fbi) = sum(prVarRcvr .* (haveSignal), 2) ...
%                                     ./ nSignals.^2;
% 
%             end
            

            % now limit to best option on each frequency band!!
            [prMeas, prVar, freqs] = deal(NaN(size(obsData.code, 1), 2));
            for s = 1:size(obsData.code, 1)
                % choose best signal on fBand 1
                sigsSF = find(obsData.fBand(s, :) < 10 & isfinite(obsData.code(s, :)));
                if ~isempty(sigsSF)
                    [prVar(s, 1), iBest] = min(prVarRcvr(s, sigsSF));
                    prMeas(s, 1) = prSmoothed(s, sigsSF(iBest));
                    freqs(s, 1) = obsData.freq(s, sigsSF(iBest));
                end
                % choose best dual frequency signal
                sigsDF = find(obsData.fBand(s, :) > 10);
                if ~isempty(sigsDF)
                    [prVar(s, 2), iBest] = min(prVarRcvr(s, sigsDF));
                    prMeas(s, 2) = prSmoothed(s, sigsDF(iBest));
                    freqs(s, 2) = obsData.freq(s, sigsDF(iBest));
                end                
            end

            % all measurement inputs are now N x 3 matricies [f_1 f_2 f_IF]

            % clean up: remove satellites without measurements
            noMeas = ~any(isfinite(prMeas) & isfinite(prVar), 2);
            prMeas(noMeas, :) = [];
            prVar(noMeas, :) = [];
            satIds(noMeas) = [];
            freqs(noMeas, :) = [];
            
            % Step 1c: limit to 1 single frequency (SF) measurement
            % NOTE: that's a bad idea. TGD, iono corrections expect L1
%             chooseSecondSF = isfinite(prSmoothed(:, 2)) ...
%                            & (prVarRcvr(:, 2) < prVarRcvr(:, 1) ...
%                               | isnan(prSmoothed(:, 1)));
%             prSF = prSmoothed(:, 1);
% %             prSF(chooseSecondSF) = prSmoothed(chooseSecondSF, 2);
%             prVarSF = prVarRcvr(:, 1);
% %             prVarSF(chooseSecondSF) = prVarRcvr(chooseSecondSF, 2);
%             
%             prMeas = [prSF, prSmoothed(:, end)];
%             prVar = [prVarSF, prVarRcvr(:, end)];
            
        end

        function ionoCorrCoeffs = getIonoCoeffs(obj, ephEp)
            %Extracts the iono correction coefficients for the right day.
            %   
            %   ionoCorrCoeffs = obj.getIonoCoeffs(epoch)
            %
            %   For a given epoch extracts the GPS iono correction
            %   coefficients for the right day.

            % pull out broadcast ephemeris struct
            eph = obj.satEph.BEph;

            [doy, yr] = navsu.time.jd2doy(navsu.time.epochs2jd(ephEp));
            ephDay = eph.doy == floor(doy) & eph.year == yr;

            % do we have GPS correction parameters for this day?
            if sum(ephDay) ~= 1 ...
                || isempty(eph.gps) ...
                || find(ephDay) > size(eph.gps.ionoCorrCoeffs, 1)
                % don't have ephemeris for this day
                ionoCorrCoeffs = NaN(1, 8);
            else
                ionoCorrCoeffs = eph.gps.ionoCorrCoeffs(ephDay, :);
            end
            
        end
        
        function [codeObservables, fn] = findRnxCodeObs(~, rnxStruct, sats, ep)
            % Returns a logical matrix indicating which observables in the
            % rinex struct contain code observables of each frequency.
            % Used by obj.readRinexData. Needs input struct in the same
            % format.
            %   [codeObservables, fn] = rnxCodeObs(~, rnxStruct, sats, ep)

            % identify code measurements among observables
            fn = fieldnames(rnxStruct.meas);
            codeMeas = cellfun(@(x) strcmp(x(1), 'C') && length(x) == 3, fn);

            % identify which observable has data for which constellation
            constRows = rnxStruct.constInds(:) == unique(rnxStruct.constInds(:))';
            nConsts = size(constRows, 2);
            hasConstMeas = false(length(fn), nConsts);
            for c = 1:nConsts
                cRows = intersect(find(constRows(:, c)), sats);
                hasConstMeas(:, c) = cellfun(@(x) ...
                    any(isfinite(rnxStruct.meas.(x)(cRows, ep))), fn);
            end
            % how many signals per constellation
            codeObservables = hasConstMeas & codeMeas;

        end

        function [obsData, SigRNM] = buildIFobs(obj, satIds, obsData)
            % Calculate the iono-free (IF) combinations of the measurements
            % contained in the obsData struct.
            % 
            % obsIF = obj.buildIFobs(obsData)
            % 
            % Inputs:
            %   obsData   struct of obs data. See obj.preprocessMeas(...)
            
            % count number of iono-free combinations
            nSigs = size(obsData.code, 2);

            % determine possible iono-free signal combinations 
            ifCombs = false(nSigs-1, nSigs);
            for sI = 1:nSigs-1

                % for valid IF combination:
                %   - both fBands finite: isfinite(obsData.fBand(:, sI)) &
                %                         isfinite(obsData.fBand(:, sI+1:end))
                %   - fBands not equal: obsData.fBand(:, sI) ~=
                %                       obsData.fBand(:, sI+1:end)

                ifCombs(sI, sI+1:nSigs) = ...
                    any(isfinite(obsData.fBand(:, sI)) ...
                      & isfinite(obsData.fBand(:, sI+1:end)) ...
                      & obsData.fBand(:, sI) ~= obsData.fBand(:, sI+1:end), 1);
            end

            % get recever tracking noise + multipath of code, carrier meas
            SigRNM.code = obj.codeMeasSigma(satIds, obsData.freq).^2;
            SigRNM.carrier = obj.carrierMeasSigma(obsData.freq).^2;

            % build IF combinations
            for sI = find(any(ifCombs, 2))'
                gamma = (obsData.freq(:, sI) ./ obsData.freq(:, ifCombs(sI, :))).^2;
                % build iono-free code and carrier measurements
                for o = {'code', 'carrier'}
                    oIF = obj.buildIFmeas(obsData.(o{1})(:, sI), ...
                                          obsData.(o{1})(:, ifCombs(sI, :)), ...
                                          gamma);
                    obsData.(o{1}) = [obsData.(o{1}), oIF];
                    % set receiver noise + multipath variances as well
                    if nargout > 1
                        vIF = obj.buildIFvar(SigRNM.(o{1})(:, sI), ...
                                             SigRNM.(o{1})(:, ifCombs(sI, :)), ...
                                             gamma);
                        SigRNM.(o{1}) = [SigRNM.(o{1}), vIF];
                    end
                end

                % determine carrier lock time of IF combinations
                obsData.tLock = [obsData.tLock, ...
                                 min(obsData.tLock(:, sI), obsData.tLock(:, ifCombs(sI, :)))];

                % create IF signal identifiers
                for IFi = find(ifCombs(sI, :))
                    newRnx = repmat({''}, size(obsData.rnxCode(:, sI)));
                    haveDF = ~any(cellfun(@isempty, obsData.rnxCode(:, [sI, IFi])), 2);
                    joinedCodes = join([repmat({'IF'}, size(newRnx)), ...
                                        obsData.rnxCode(:, [sI, IFi])], '');
                    newRnx(haveDF) = joinedCodes(haveDF);
                    obsData.rnxCode = [obsData.rnxCode, newRnx];
                end

                % make up IF frequency bands
                for IFi = find(ifCombs(sI, :))
                    newBand = obsData.fBand(:, [sI, IFi]) * [10; 1];
                    obsData.fBand = [obsData.fBand, newBand];
                end

                % fill remaining properties (CN0, Doppler) with NaNs
                for fn = fieldnames(obsData)'
                    if size(obsData.(fn{1}), 2) == nSigs
                        % no iono-free values were added yet
                        obsData.(fn{1}) = [obsData.(fn{1}), NaN(size(oIF))];
                    end
                end

                
            end
        end

        function ifMeas = buildIFmeas(~, measF1, measF2, gamma)
            % Builds the iono-free measurement combination for two
            % measurements (code or carrier) and the ratio of their
            % frequencies, gamma.
            % 
            % ifMeas = obj.buildIFmeas(measF1, measF2, gamma)
            %   
            % where gamma = (F1/F2).^2

            ifMeas = (gamma.*measF1 - measF2) ./ (gamma-1);
        end
        
        function varIF = buildIFvar(~, measVarF1, measVarF2, gamma)
            % Calculate the variance of the iono-free (IF) combination of 
            % measurements with given variances and the ratio of their
            % frequencies, gamma.
            % 
            % varIF = obj.buildIFobs(measVarF1, measVarF2, gamma)
            % 
            % measVar1  N x 1 matrix of meas variances (code or carrier) in
            %           meter^2 on frequency 1
            % measVar2  N x M matrix of meas variances (code or carrier) in
            %           meter^2 on frequency 2
            % gamma     N x M matrix of squared frequency ratios (F1/F2).^2

            varIF = (gamma.^2.*measVarF1 + measVarF2) ./ (gamma-1).^2;
        end

        function [prSmoothed, prVarRcvr] = carrierSmoothing(obj, satIds, epoch, obsData, SigRNM)
            % Performs carrier smoothing on the observables.
            %   Returns a matrix of smoothed range measurements based on
            %   passed code and carrier values. Also returns a matrix of
            %   the approximate receiver noise on the smoothed signals.
            %   
            %   [prSmoothed, prVarRcvr] = obj.carrierSmoothing( ...
            %                              satIds, epoch, obsData, SigRNM)

            % Need one carrier smoother per signal, per constellation
            % Each smoother is identified by constelleation and rnx sig id


            % initialize smoothed values as equal to code only
            prSmoothed = obsData.code;
            prVarRcvr = SigRNM.code;

            if obj.useCarrierSmoothing
                % find all the received signals
                haveSigs = cellfun(@(x) ~isempty(x), obsData.rnxCode) ...
                         & isfinite(obsData.code);

                for c = unique(obj.satConstId(satIds))'

                    % from which satellites did we receive signals?
                    constSigs = obj.satConstId(satIds) == c & haveSigs;
                    uniqueSigs = unique(obsData.rnxCode(constSigs));

                    for s_i = 1:length(uniqueSigs)
                        % call the right smoother for each constellation,
                        % each signal type
                        
                        % grab the signal name
                        sigName = uniqueSigs{s_i};

                        % get indices of these signals
                        sigIds = constSigs ...
                               & strcmp(sigName, obsData.rnxCode);

                        if any(sum(sigIds, 2) > 1)
                            warning('I found multiple %s "%s" signals! Skipping.', ...
                                navsu.svprn.convertConstIndName(c), sigName);
                            warning('Check Rinex signal identifyers.');
                            continue
                        end
                        
                        smootherId = strcmp(sigName, {obj.CS.signal}) ...
                                   & c == [obj.CS.const];
    
                        if ~any(smootherId)
                            % need to add a smoother for this signal
                            if startsWith(sigName, 'IF')
                                % set large smoothing time for iono-free combo
                                tMax = obj.smoothingConstantIF;
                            else
                                tMax = obj.smoothingConstant;
                            end
                            obj.CS(end+1) = navsu.lsNav.CarrierSmoother(sigName, ...
                                                                        c, ...
                                                                        obj.numSats, ...
                                                                        tMax);
                            % use this new smoother
                            smootherId = [smootherId, true]; %#ok
                        end
                        [prSmoothed(sigIds), prVarRcvr(sigIds)] = ...
                            obj.CS(smootherId).smoothen( ...
                            satIds(any(sigIds, 2)), epoch, ...
                            structfun(@(x) x(sigIds), obsData, 'UniformOutput', false), ...
                            structfun(@(x) x(sigIds), SigRNM, 'UniformOutput', false));
                    end
    
                end
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
            
            if obj.usePreciseProducts
                % Propagate the precise ephemeris
                % need: satPos, satVel, satClk.bias, satClk.drift, posVar,
                % TGD
                obj.updateSatDataFromPProd(satIds(n2p), epoch(n2p));
            else
                % Propagate the broadcast.
                obj.updateSatDataFromBrdc(satIds(n2p), epoch(n2p));
            end
            
        end
        
        function updateSatDataFromBrdc(obj, s2update, epoch)
            %Propagate the navigation broadcast for selected satellites to
            %chosen epoch. Update internal memory accordingly.
            
            % here I could check the memory to see what I need to propagate                
                
            % propagate navigation broadcast to desired epoch
            satInfo = navsu.geo.propNavMsg(obj.satEph.BEph, ...
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
            %   s2update    index of satellites to compute data for
            %   epoch       epoch at which to compute the data for
            
            N = sum(s2update>0);
            [obj.satPos(s2update, :), obj.satVel(s2update, :)] = ...
                obj.satEph.propagate(obj.satPRN(s2update), ...
                                     obj.satConstId(s2update), ...
                                     epoch.*ones(N, 1));
            svClock = obj.satEph.clock(obj.satPRN(s2update), ...
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
            constIds = obj.getActiveConstId(obj.satConstId(satIds));
            
            ttx = sum([measEpoch*ones(sum(satIds > 0), 1), ...
                - obj.tBias(constIds) / navsu.constants.c, ...
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
            
            % check measurement validity
            measExclude = obj.satEl(satIds) < obj.elevMask ...
                        | ~isfinite(prhat) ...
                        | ~isfinite(measVar);
            satIds(measExclude) = [];
            prhat(measExclude) = [];
            measVar(measExclude) = [];
            
            % check again which states to really compute
            usedConsts = unique(obj.satConstId(satIds))';
            outStates = [1:3, 3+find(ismember(usedConsts, outConsts))];
            
            % initialize outputs
            N = length(measExclude);
            dx = NaN(nOutStates, 1);
            Rout = NaN(nOutStates);
            prr = NaN(N, 1);
            P = NaN(N);
            DOP = NaN(nOutStates, nOutStates);
            
            % make sure there's enough satellites remaining
            if length(outStates) > length(prhat)
                return
            end

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
                prr(~measExclude) = prhat - G*dx(outStates);
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
        
        function obsData = frequencyMask(~, obsData, fBands)
            % Limit obsData to signals on specific frequency bands
            % obsData = obj.frequencyMask(obsData, fBands)

            % choose frequency bands
            if max(fBands) <= size(obsData.code, 2)
                obsData = structfun(@(x) x(:, fBands), obsData, ...
                                    'UniformOutput', false);
            else
                warning('Illegal frequency band selection.');
            end

        end

        function usedConstIdx = getActiveConstId(obj, constIds)
            % Compute the index among the used constellations for a const
            % id.

            if isrow(constIds)
                % usually comes as a column vector
                constIds = constIds';
            end

            if ~all(constIds, obj.activeConsts)
                % check if only legal constellations passed
                error('Did not expect these constellations!');
            end
            constMatches = find(obj.activeConsts == constIds');
            % compute index among used consts for each const id
            usedConstIdx = constMatches(:) ...
                         - obj.numConsts*(0:length(constIds)-1)';
        end

        function updateClkBias(obj, measEpoch)
            % Updates the clock bias estimate at the given epoch based on
            % previous estimates of clock bias and clock bias rate.
            consts = isfinite(obj.tBiasRate);
            if isfinite(obj.tPosSol)
                obj.tBias(consts) = obj.tBias(consts) ...
                      + (measEpoch - obj.tPosSol) .* obj.tBiasRate(consts);
            end
        end
    end
end

