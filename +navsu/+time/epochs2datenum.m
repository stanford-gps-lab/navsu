function tDatenum = epochs2datenum(epochs)
% epochs2datenum
% DESCRIPTION:
%   Convert from GPS epochs(seconds since start of GPS time) to MATLAB
%   datenum
% INPUT:
%   epochs = Nx1 vector of GPS epochs
%
% OUTPUT:
%   datenums = Nx1 vector of MATLAB datenums
%
% See also: navsu.time.datenum2epochs
tDatenum = datenum(navsu.time.epochs2cal(epochs,1));


end