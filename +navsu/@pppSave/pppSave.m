classdef pppSave < handle

    properties
        posSave
        velSave
        attSave
        attEnuSave
        biasSave
        epochSave
        clockSave
        INDS_STATE
        covSave
        measType
        
        gnssData = struct('range',[],'doppler',[],'epochs',[],'ind',[],'sig',[],...
            'freqs',[],'PRN',[],'constInds',[],'rnxCode',[],'el',[]);
        
        figCov
        figResidsSummary
        figResids
        figResidsSummaryGlo
        figResidsEl 
        
        tdxSmall
        stateSaveSmall
        covSaveSmall
        epochsSmall
        covEnuSmall
        plSmall
        plLocSmall
        
        tdxFull
        stateSaveFull
        covSaveFull
        epochsFull
        tdx2
        covEnuFull
        plFull
        plLocFull
        
        % Stuff for removed measurements
        figRemoved
        epochsRemoved 
        measRemoved
        indRemoved
        
        % Solution separation info if necessary
        
    end


    methods
        function obj = pppSave(epochsFull,nConst,obsGnss,epochsGnss,...
                INDS_STATE,measType,varargin)
            
            p = inputParser;
            
            p.addParameter('gnssEpochsOnly',true);
            
            % parse the results
            parse(p, varargin{:});
            res        = p.Results;
            gnssEpochsOnly  = res.gnssEpochsOnly;  % Index of full state saving matrices
            
            if nargin < 1
                % if no values are provided (this is the only option right
                % now)
                obj.posSave    = [];
                obj.velSave    = [];
                obj.attSave    = [];
                obj.attEnuSave = [];
                obj.biasSave   = [];
                obj.epochSave  = [];
                obj.clockSave  = [];
                obj.covSave    = [];
                obj.gnssData.rangeResids = [];
                obj.gnssData.epochs      = [];
                obj.INDS_STATE = [];
                
            else
                if gnssEpochsOnly
                    epochsFull = epochsGnss;
                end
                nEpochsFull = length(epochsFull);
                
                % if values are provided
                obj.posSave = nan(nEpochsFull,3);
                obj.velSave = nan(nEpochsFull,3);
                obj.attSave = nan(nEpochsFull,3);
                obj.attEnuSave = nan(nEpochsFull,3);
                obj.biasSave = nan(nEpochsFull,6);
                obj.epochSave = epochsFull;
                obj.measType  = measType;
                obj.clockSave = nan(nEpochsFull,nConst);
                obj.INDS_STATE           = INDS_STATE;
                obj.covSave   = nan(nEpochsFull,INDS_STATE.FLEX_STATE_MIN-1);
                
                
                % Range information to be saved
                rangeMeas = obsGnss.range;
                
                obj.gnssData.range.resids      = nan(size(rangeMeas.obs));
                
                obj.gnssData.range.ind         = rangeMeas.ind;
                obj.gnssData.range.sig         = rangeMeas.sig;
                obj.gnssData.range.freqs       = rangeMeas.freqs;
                obj.gnssData.range.PRN         = rangeMeas.PRN;
                obj.gnssData.range.constInds   = rangeMeas.constInds;
                obj.gnssData.range.rnxCode     = rangeMeas.rnxCode;
                
                
                % General satellite info
                obj.gnssData.PRN       = obsGnss.PRN;
                obj.gnssData.constInds = obsGnss.constInds;
                obj.gnssData.freqs     = obsGnss.freqs;
                obj.gnssData.el        = nan(size(obj.gnssData.range.resids,2),size(obj.gnssData.range.resids,3));
                obj.gnssData.epochs    = obsGnss.epochs;
                
                  % Doppler information to be saved
                dopplerMeas = obsGnss.doppler;
                
                obj.gnssData.doppler.resids      = nan(size(dopplerMeas.obs));
                
                obj.gnssData.doppler.sig         = dopplerMeas.sig;
                obj.gnssData.doppler.freqs       = dopplerMeas.freqs;
                obj.gnssData.doppler.PRN         = dopplerMeas.PRN;
                obj.gnssData.doppler.constInds   = dopplerMeas.constInds;
                obj.gnssData.doppler.rnxCode     = dopplerMeas.rnxCode;
                
                obj.figCov = figure('visible','off');
                obj.figResids = figure('visible','off');
                obj.figResidsSummary = figure('visible','off');
                obj.figResidsSummaryGlo = figure('visible','off');
                obj.figRemoved = figure('visible','off');
                obj.figResidsEl = figure('visible','off');
                
                obj.tdxSmall = [];
                obj.stateSaveSmall = [];
                obj.covSaveSmall = [];
                obj.epochsSmall = [];
                obj.epochsFull = epochsFull;
                obj.stateSaveFull = nan(19,length(epochsFull));
                obj.covSaveFull   = nan(INDS_STATE.FLEX_STATE_MIN-1,nEpochsFull);
                obj.covEnuFull    = nan(3,nEpochsFull);
                obj.plFull = nan(3,nEpochsFull);
                obj.plLocFull = nan(3,nEpochsFull);
                
                obj.epochsRemoved = [];
                obj.measRemoved   = [];
                obj.indRemoved    = [];
            end
        end
    end


    % function signatures
    methods
        saveNewValues(obj,filter,PARAMS,varargin)
         [rangeResids, doppResids, elFull, azFull]  = saveResids(obj,measMat,residsPost,epoch,el,az,prnConstInds)
        plotResidSummary(obj,varargin)
        [epochs,posTrue,posEst,pl,errEnu,plLoc] = plotAgainstTruth(obj,filenameTruth,PARAMS,varargin)
        midRunPlot(obj)
        plotResids(obj)
        closeSaveState(obj)
        [measRemove,epoch] = saveMeasRemoved(obj,epoch,measLow,measResids,measSlip)
        plotRemoved(obj)
        plotSol(obj,varargin)
    end


end
