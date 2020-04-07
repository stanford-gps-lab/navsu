function [freqs, freqInds] = mapSignalFreq(freqDes,prns,constInds,jd)
% mapSignalFreq
% DESCRIPTION:
%   Map from the signal number used in RINEX observation codes (1, 2, 5, 6,
%   7, 8)  and the constellation to the actual frequency of the signal
% INPUTS:
%   freqDes   - a 5xN matrix, where N is the number of desired signals output,
%               and 5 is the number of constellations we consider (GRECS)
%   prns      - vector length M, list of PRNs
%   constInds - vector length M, list of indices indicating constellation (GRECS) of PRNs
%   jd        - julian date, if necessary, to indicate GLONASS freq
%
% OUTPUTS:
% freqs:     MxN matrix of frequencies associated with each signal 
%
% See also: navsu.svprn.prn2svn, navsu.svprn.prn2svn


fRinexLabels = [1 2 5 6 7 8]';

fTable = [1575.42	1602	1575.42     nan      nan	1575.42;
          1227.6	1246	nan         1561.098 nan	1227.6;
          1176.45	1201	1176.45     nan      nan    1176.45;
          nan	    nan     1278.75     1268.52	 nan	1278.75;
          nan       nan     1207.14     1207.14	 nan	nan;
          nan       nan     1191.795	nan      nan	nan];

% GLONASS FDMA channel separation
gloChannels = [0.5625 0.4375 0.4375];      

freqs = nan(length(prns),size(freqDes,2));
freqInds = nan(length(prns),size(freqDes,2));
for pdx = 1:length(prns)
    consti = constInds(pdx);
    
    for idx = 1:size(freqDes,2)
        fDesi = freqDes(pdx,idx);
        
        if isnan(fDesi)
            continue;
        end
        
        freqi = fTable(find(fRinexLabels == fDesi),consti);
        
        if consti == 2 && fDesi <= 3
            % NEED TO ACCOUNT FOR GLONASS FDMA OFFSET
            freqChani = utility.svprn.prn2FreqChanGlonass(prns(pdx),jd);
            freqi = freqi+gloChannels(fDesi)*freqChani;
        end
        freqs(pdx,idx) = freqi;
        
        freqInds(pdx,idx) = freqDes(pdx,idx);
    end
end

freqs = freqs*1e6;
end




















