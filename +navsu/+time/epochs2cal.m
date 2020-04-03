function [yr,mn,dy,hr,min,sec] = epochs2cal(epochs,vectorFlag)

if nargin < 2
    vectorFlag = 0;
end


% Ensure correct direction
if size(epochs,2) == 1
    epochs = epochs';
end

% Convert from GPS epochs to calendar date and time
[weeks,tows] = navsu.time.epochs2gps(epochs);
towDays = 86400*floor(tows/86400);

jds = navsu.time.gps2jd(weeks,towDays);
[yr,mn,dy] = navsu.time.jd2cal(jds);

hr = floor((tows-towDays)/3600);
min = floor((tows-towDays-hr*3600)/60);
sec = ((tows-towDays-hr*3600-min*60));

% Output as single vector
if vectorFlag
    yr = [yr' mn' dy' hr' min' sec'];
    mn = [];
    dy = [];
    hr = [];
    min = [];
    sec = [];
end

end