function [snrMaskRange, snrMaskDoppler] = snrMask(obs,threshs)
% snrMask  given a GNSS observation and thresholds for each frequency,
% produce masking matrices for the range and doppler measurements
% (obs.range.obs and obs.doppler.obs).  See also
% navsu.estimators.pppFilter.handleGnssMeas
%
% INPUTS: 
%   threshs  [1x3] vector of Cn0 thresholds for signals 1, 2, and 3
%   obs      Single GNSS obs
% OUTPUTS: 
%   snrMaskRange logical matrix, true = SNR check passed for obs.range.obs
%                elements
%   snrMaskDoppler logical matrix, true = SNR check passed for obs.doppler.obs
%                  elements
% 

% Range mask of SNR
sigs = unique(obs.range.sig);
snrRange = zeros(size(obs.range.obs,1),size(obs.range.obs,2),2);
snrThresh = zeros(size(obs.range.obs,1),size(obs.range.obs,2),2);
for idx =1:length(sigs)
    sigsi = [mod(sigs(idx),100) floor(sigs(idx)/100)];
    
    for jdx = 1:2

        sigi = sigsi(jdx);
        
        rowRange = find(obs.range.sig(:,1) == sigs(idx));
        if sigi == 0
            snrRange(rowRange,:,jdx) = Inf;
            
            continue;
        end
        
        rowSnr = find(obs.snr.sig(:,1) == sigi);

        snrRange(rowRange,:,jdx) = repmat(obs.snr.obs(rowSnr,:),length(rowRange),1);
        snrThresh(rowRange,:,jdx) = threshs(sigi);
    end
end

snrMaskBase = snrRange > snrThresh | repmat(isnan(obs.range.freqs),1,1,2);
snrMaskRange = all(snrMaskBase,3);



% Doppler mask of SNR

sigs = unique(obs.doppler.sig);
snrDoppler = zeros(size(obs.doppler.obs,1),size(obs.doppler.obs,2),2);
snrThresh = zeros(size(obs.doppler.obs,1),size(obs.doppler.obs,2),2);
for idx =1:length(sigs)
    sigsi = [mod(sigs(idx),100) floor(sigs(idx)/100)];
    
    for jdx = 1:2
        
        sigi = sigsi(jdx);
        
        rowRange = find(obs.doppler.sig(:,1) == sigs(idx));
        if sigi == 0
            snrDoppler(rowRange,:,jdx) = Inf;
            
            continue;
        end
        
        rowSnr = find(obs.snr.sig(:,1) == sigi);

        snrDoppler(rowRange,:,jdx) = repmat(obs.snr.obs(rowSnr,:),length(rowRange),1);
        snrThresh(rowRange,:,jdx) = threshs(sigi);
    end
end

snrMaskBase = snrDoppler > snrThresh | repmat(isnan(obs.doppler.freqs),1,1,2);
snrMaskDoppler = all(snrMaskBase,3);


end