function [freqs, freqDes] = mapSignalFreq(freqDes, prns, constInds, jd)
% [freqs, freqDes] = mapSignalFreq(freqDes, prns, constInds, jd)
% DESCRIPTION:
%   Map from the signal number used in RINEX observation codes (1, 2, 5, 6,
%   7, 8)  and the constellation to the actual frequency of the signal
% INPUTS:
%   freqDes   - a MxN matrix of Rinex frequency band designators, where N 
%               is the number of desired signals output, and M is the
%               number of PRNs we consider
%   prns      - vector length M, list of PRNs
%   constInds - vector length M, list of indices indicating constellation (GRECS) of PRNs
%   jd        - julian date, if necessary, to indicate GLONASS freq
%
% OUTPUTS:
% freqs:     MxN matrix of frequencies associated with each signal
%
% See also: navsu.svprn.prn2svn, navsu.svprn.prn2svn

%% parse inputs

% get number of satellites/frequencies
M = max([size(freqDes, 1), length(prns), length(constInds)]);

if size(freqDes, 1) == 1
    freqDes = repmat(freqDes, M, 1);
end

% want prns as M x 1 vector
if length(prns) == 1
    prns = prns * ones(M, 1);
elseif isrow(prns)
    prns = prns';
end
% same for constInds
if length(constInds) == 1
    constInds = constInds * ones(M, 1);
elseif isrow(constInds)
    constInds = constInds';
end

fRinexLabels = [1 2 3 5 6 7 8]';

% to be accessed by (fDes == fRinexLabels, const)
%         GPS       GLONASS  GAL         BDS        QZSS    SBAS
fTable = [1575.42	1602	 1575.42     1575.42    1575.42 1575.42;  % 1
          1227.6	1246	 nan         1561.098   1227.6  nan;      % 2
          nan       1202.025 nan         nan        nan     nan; % 3
          1176.45	nan      1176.45     1176.45    1176.45 1176.45;  % 5
          nan       1248.06  1278.75     1268.52    1278.75	nan;      % 6
          nan       nan      1207.14     1207.14    nan	    nan;      % 7
          nan       nan      1191.795	 1191.795   nan	    nan];     % 8

% GLONASS FDMA channel separation
gloChannels = [0.5625 0.4375];

% find index of each frequency label
logFreqLabelsIdx = fRinexLabels == freqDes(:)';
labelExists = any(logFreqLabelsIdx', 2);
freqLabelsIdx = find(logFreqLabelsIdx) - (find(labelExists)-1)*length(fRinexLabels);
% now get linear indices for fTable
correspondingConsts = repmat(constInds, size(freqDes, 2), 1);
linFreqIdx = sub2ind(size(fTable), freqLabelsIdx, correspondingConsts(labelExists));

% now retrieve frequency for each signal where I have one
freqs = nan(size(freqDes));
freqs(labelExists) = fTable(linFreqIdx);

% NEED TO ACCOUNT FOR GLONASS FDMA OFFSET
isGloSat = constInds == 2;
hasGloChannels = isGloSat & freqDes < 3 & freqDes > 0;

if any(hasGloChannels, 'all')
    if ~exist('jd', 'var')
        error('Need a julian date to get GLONASS freq. channel!')
    end
    % get each satellite's frequency channel
    freqChani = zeros(size(constInds));
    freqChani(isGloSat) = navsu.svprn.prn2FreqChanGlonass(prns(isGloSat), jd);
    % and the offset of each quieried frequency
    channelOffsets = zeros(size(freqDes));
    channelOffsets(hasGloChannels) = gloChannels(freqDes(hasGloChannels));

    % add offset to center frequencies
    freqs = freqs + channelOffsets.*freqChani;
end

% convert to MHz
freqs = freqs*1e6;
end

