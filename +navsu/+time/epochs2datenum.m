function tDatetime = epochs2datenum(epochs)
% epochs2datenum
% DESCRIPTION:
%   Convert from GPS epochs(seconds since start of GPS time) to MATLAB
%   datetime
% INPUT:
%   epochs = Nx1 vector of GPS epochs
%
% OUTPUT:
%   datenums = Nx1 vector of MATLAB datetime
%
% See also: navsu.time.datenum2epochs
tDatetime = datetime(navsu.time.epochs2cal(epochs,1));

end
