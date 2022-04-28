function [dtime] = epochs2datetime(epochs)
% epochs2datetime
% DESCRIPTION:
%   Convert from GPS epoch (seconds since start of GPS time) to matlab
%   datetime format
% INPUT:
%   epochs = Nx1 list of GPS epochs 
% OUTPUT:
%   dtime  =  [Nx1] list of datetimes
%
% See also: navsu.time.epochs2cal


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

if ~isempty(yr)
    for idx = size(yr,1):-1:1
        dtime(:,idx) = datetime([yr(idx,:)' mn(idx,:)' dy(idx,:)' hr(idx,:)' min(idx,:)' sec(idx,:)']);
    end
else
    dtime = [];
end
end