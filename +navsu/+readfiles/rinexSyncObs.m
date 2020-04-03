function [time, date, obsStruc, max_int] = ...
          rinexSyncObs(time_i, week_i, date_i, obsStruc_i, interval)

% SYNTAX:
%   [time_ref, time, week, date, pr1, ph1, pr2, ph2, dop1, dop2, snr1, snr2] = ...
%   sync_obs(time_i, week_i, date_i, pr1_i, ph1_i, pr2_i, ph2_i, dop1_i, dop2_i, snr1_i, snr2_i, interval);
%
% INPUT:
%   time_i = receiver seconds-of-week
%   week_i = GPS week
%   date_i = date (year,month,day,hour,minute,second)
%   pr1_i = code observation (L1 carrier)
%   ph1_i = phase observation (L1 carrier)
%   pr2_i = code observation (L2 carrier)
%   ph2_i = phase observation (L2 carrier)
%   dop1_i = Doppler observation (L1 carrier)
%   dop2_i = Doppler observation (L2 carrier)
%   snr1_i = signal-to-noise ratio (L1 carrier)
%   snr2_i = signal-to-noise ratio (L2 carrier)
%   interval = observation time interval [s]
%
% OUTPUT:
%   time_ref = reference seconds-of-week
%   time = receiver seconds-of-week
%   week = GPS week
%   date = date (year,month,day,hour,minute,second)
%   pr1 = code observation (L1 carrier)
%   ph1 = phase observation (L1 carrier)
%   pr2 = code observation (L2 carrier)
%   ph2 = phase observation (L2 carrier)
%   dop1 = Doppler observation (L1 carrier)
%   dop2 = Doppler observation (L2 carrier)
%   snr1 = signal-to-noise ratio (L1 carrier)
%   snr2 = signal-to-noise ratio (L2 carrier)
%
% DESCRIPTION:
%   Synchronize different sets of observations. Zeros where not available.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

obsFields = fieldnames(obsStruc_i);

%number of satellite slots
nSatTot = size(obsStruc_i.(obsFields{1}),1);

%number of observation datasets (e.g. number of read RINEX files)
nObsSet = size(obsStruc_i.(obsFields{1}),3);

%find min and max time tags (in common among all observation datasets)
time_i_nan = time_i;

tNoMeas = true(size(obsStruc_i.(obsFields{1}),2),1);
for idx =1 :length(obsFields)
    tNoMeas = tNoMeas & ~(permute(sum(obsStruc_i.(obsFields{idx}),1),[2 1 3])==0);
end
time_i_nan(tNoMeas) = NaN; %set NaN to epochs which don't have any pseudorange
min_time_prog = min(min(time_i_nan,[],1));

%find the largest interval
max_int = max(interval(:));
%max_int = 30;

%define the reference time
time_ref = unique(time_i);
time_ref(isnan(time_ref)) = [];

tow_ref = mod(time_ref,60*60*24*7);
tow_ref=roundmod(tow_ref,max_int);

%number of reference epochs
ref_len = length(tow_ref);

%create containers
time = zeros(ref_len, 1, nObsSet);
week = zeros(ref_len, 1, nObsSet);
date = zeros(ref_len, 6, nObsSet);

for idx = 1:length(obsFields)
    obsStruc.(obsFields{idx}) = zeros(nSatTot, ref_len, nObsSet);
end

time_prog = time_i - min_time_prog; % substract the first element to reduce the magnitude of all the values
time_ref_prog = time_ref - min_time_prog;

for s = 1 : nObsSet
    [~, idx_t, idx_z] = intersect(roundmod(time_ref_prog, max_int), roundmod(time_prog(:,1,s), interval(s)));    

    
    time(:,s) = time_ref;
    week(:,s) = navsu.time.epochs2gps(time_ref);
    date(:,:,s) = navsu.time.epochs2cal(time_ref,1);
    
    for idx = 1:length(obsFields)
        obsStruc.(obsFields{idx})(:,idx_t,s) = obsStruc_i.(obsFields{idx})(:,idx_z,s);
    end   
    
end

end

function rx = roundmod(x,y)

rx = round(x./y).*y;

end

