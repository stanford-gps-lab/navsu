function [Pxyzt, Vxyzt] = interpolateEphemeris(prn,Peph,Clck,epochDayStart, epochDayEnd,settings)

% Pick off desired PRN
Pposition    = Peph.position(Peph.prns == prn,:);
Pvelocity    = Peph.velocity(Peph.prns == prn,:);
Pclock_bias  = Peph.clock_bias(Peph.prns == prn,:);
Pclock_drift = Peph.clock_drift(Peph.prns == prn,:);
Pepochs      = Peph.epochs(Peph.prns == prn,:);
Pevent       = Peph.event(Peph.prns == prn,:);
Cepochs      = Clck.Cepochs;
Cclock       = Clck.Cclk(prn,:);
Cclock_sig   = Clck.Cclk_sig(prn,:);

% Set up interpolation constants
% n     = settings.nPolyFit; % number of nominal points to use for polynominal fit
n     = 9;                  % RW - 8/22 test, there should be 9 points fed into polyinterp
pfit  = settings.pfit;     % order of fit for position interpolation
cdfit = settings.cdfit;    % order of fit for clock interpolation

process_interval = Cepochs(2) - Cepochs(1);
n_epochs_day = 86400/process_interval;
cdecimate = Peph.Epoch_interval/process_interval; %how many high rate epochs per low rate epoch

% rcidx = (1:n) - n/2;                % relative indices to fit
rcidx = (1:n) - floor(n/2);         % RW - 8/22 test, indices can't be noninteger
rpidx = rcidx;                      % relative indices for precise orbit file
rtf   = rcidx*Peph.Epoch_interval;  % relative time steps for fits
rcidx = rcidx*cdecimate;            % relative indices for high rate clock data
rcjdx = (1:(cdecimate-1));          % relative indices for high rate clock data
rpjdx = rcjdx;                      % relative indices for high rate interpolated data
rtjdx = rcjdx'*process_interval;    % relative time steps for high rate interpolated data

chi2 = NaN(n_epochs_day*Peph.NumSV/cdecimate, 4);
pvar = ones(size(rpidx))'*1e-7;     % supposed to scale the chi^2 statistic output?

% Initialize
Pxyzt = nan(n_epochs_day,4);
Vxyzt = nan(n_epochs_day,4);

% Find start and end indices
idx_start = find(Pepochs >= epochDayStart); 
idx_start = idx_start(1);

idx_end = find(Pepochs < epochDayEnd); 
idx_end = idx_end(end);

for idx = idx_start:idx_end
    Pxyzti = NaN(cdecimate, 4);
    Vxyzti = NaN(cdecimate, 4);
    
    % set first value to precise file values
    Pxyzti(1,1:3) = Pposition(idx,:);
    Pxyzti(1,4) = Pclock_bias(idx);
    Vxyzti(1,1:3) = Pvelocity(idx,:);
    Vxyzti(1,4) = Pclock_drift(idx);
    
    % interpolate precise orbit values
    % polyinterp called 3 times to calculate the interpolated values at a
    % higher rate of 30 seconds one coordinate at a time 
    % Pxyzti(1 + rcjdx, 1) output should be a 9x1

%     [Pxyzti(1 + rcjdx, 1), Vxyzti(1 + rcjdx, 1), chi2(idx,1), p] = ...
%         navsu.geo.polyinterp(rtf', Pposition(rpidx + idx,1), pfit, rtjdx, ...
%         Pevent(rpidx +idx), pvar);

    [Pxyzti(1 + rcjdx, 1), Vxyzti(1 + rcjdx, 1), chi2(idx,1)] = ...
        navsu.geo.polyinterp(rtf', Pposition(rpidx + idx,1), pfit, rtjdx, ...
        Pevent(rpidx +idx), pvar);

%     [Pxyzti(1 + rcjdx, 1), Vxyzti(1 + rcjdx, 1)] = ...
%         navsu.geo.pephInterp(Peph, prn, Pposition(rpidx + idx,1), pfit, rtjdx, ...
%         Pevent(rpidx +idx), pvar);  % RW - 8/7/21

    [Pxyzti(1 + rcjdx, 2), Vxyzti(1 + rcjdx, 2), chi2(idx,2)] = ...
        navsu.geo.polyinterp(rtf', Pposition(rpidx + idx,2), pfit, rtjdx, ...
        Pevent(rpidx +idx), pvar);
    [Pxyzti(1 + rcjdx, 3), Vxyzti(1 + rcjdx, 3), chi2(idx,3)] = ...
        navsu.geo.polyinterp(rtf', Pposition(rpidx + idx,3), pfit, rtjdx, ...
        Pevent(rpidx +idx), pvar);
    
    
    % adjust high rate clock so that it matches NGA timebase
    cedx = floor((idx-1))*cdecimate + 1;
    cprn = prn;
    cdx = cedx + [0 rcjdx];
    
    Pxyzti(1 + [0 rcjdx], 4) = Cclock(cdx);
    tmp2 = [0; rtjdx];
    % RW - what does the next line do?
    [~,Vxyzti(1 + [0 rcjdx], 4)] = ...
            navsu.geo.polyinterp(rtf', Cclock(cedx + rcidx)',cdfit,tmp2, ...
            (Pevent(rpidx +idx) | ...
            abs(Pclock_bias(rpidx + idx)) > 0.9999), ...
            (Cclock_sig(cedx + rcidx).^2)' + 1e-22);
    
%     
%     if sum(~isnan(Cclock( cedx + rcidx))) > cdfit && ...
%             any(~isnan(Pclock_bias(rpidx(rpidx <= 0) + idx))) && ...
%             any(~isnan(Pclock_bias(rpidx(rpidx > 0) + idx)))
%         % if the first point is missing from NGA replace it with CODE
%         tmp = rcjdx;
%         tmp2 = rtjdx;
%         if isnan(Pclock_bias(idx)) && ~isnan(Cclock( cedx))
%             tmp = [0 tmp];
%             tmp2 = [0; tmp2];
%         end
%         %interpolate difference between NGA and CODE
%         [Pxyzti(1 + tmp, 4), Vxyzti(1 + tmp, 4), chi2(idx,4)] = ...
%             polyinterp(rtf', Pclock_bias(rpidx + idx) - ...
%             Cclock(cedx + rcidx)',cdfit,tmp2, ...
%             (Pevent(rpidx +idx) | ...
%             abs(Pclock_bias(rpidx + idx)) > 0.9999), ...
%             (Cclock_sig(cedx + rcidx).^2)' + 1e-22);
%         %smooth CODE clock drift estimate
%         tc = (((rtf(1)+process_interval):process_interval:rtf(end))' + (rtf(1):process_interval:(rtf(end)-process_interval))')/2;
%         cc = diff(Cclock(cedx + (rcidx(1):rcidx(end)))')/process_interval;
%         fc = zeros(size(tc));
%         mc = nanmean(cc);
%         sc = nanstd(cc);
%         fc(abs(cc - mc) > 4.5*sc) = 1;
%         if sum(fc)
%             % if any points are removed, recalcualte mean and sigma
%             cc(abs(cc - mc) > 4.5*sc) = NaN;
%             mc = nanmean(cc);
%             sc = nanstd(cc);
%             fc(abs(cc - mc) > 4.5*sc) = 1;
%         end
%         vc = (Cclock_sig(cedx + (rcidx(1):rcidx(end))).^2)' + 1e-22;
%         [tmp3, ~, ~] = polyinterp(tc, cc ,cdfit,tmp2, fc, vc);
%         clock_diff(prn) = mean(Pxyzti(1 + tmp, 4));
%         drift_diff(prn) = mean(Vxyzti(1 + tmp, 4));
%         Pxyzti(1 + tmp, 4) = Pxyzti(1 + tmp, 4) + Cclock(cedx + tmp)';
%         Vxyzti(1 + tmp, 4) = Vxyzti(1 + tmp, 4) + tmp3;
%     end
    
    % Save off interpolated values
    Pxyzt((idx-idx_start)*cdecimate+1:(idx-idx_start+1)*cdecimate,:) = Pxyzti;
    Vxyzt((idx-idx_start)*cdecimate+1:(idx-idx_start+1)*cdecimate,:) = Vxyzti;
    
end % idx = idx_start:idx_end

end























