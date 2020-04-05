function epochs = datenum2epochs(datenums)
% datenum2epochs
% DESCRIPTION:
%   Convert from MATLAB datenum to GPS epochs (seconds since start of GPS
%   time)
% INPUT:
%   datenums = Nx1 vector of MATLAB datenums
%
% OUTPUT:
%   epochs = Nx1 vector of GPS epochs
%
% See also: navsu.time.epochs2datenum

[Y,M,D,H,MN,S] = datevec(datenums);

epochs = navsu.time.cal2epochs(Y,M,D,H,MN,S);

end
