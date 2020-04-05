function [yr,mn,dy,hr,min,sec] = epochs2cal(epochs,vectorFlag)
% epochs2cal
% DESCRIPTION:
%   Convert from GPS epoch (seconds since start of GPS time) to calendar
%   date
% INPUT:
%   epochs     = Nx1 list of GPS epochs 
%   vectorFlag = optional input to change outputs to an Nx6 single matrix
%              combining all outputs
% OUTPUT:
%   yr  = year      [Nx1] or [yr mn dy hr min sec] [Nx6] if vectorFlag = 1
%   mn  = month     [Nx1]
%   dy  = day       [Nx1]
%   hr  = hour      [Nx1]
%   min = minute    [Nx1]
%   sec = second    [Nx1]
%
% See also: navsu.time.cal2epochs


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