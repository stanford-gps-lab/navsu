function epochs = jpl2epochs(jplTime)
% jpl2epochs
% DESCRIPTION:
% Convert from the time used in JPL real time ephemeris products to GPS
% epochs (seconds since start of GPS time).  JPL time is very similar to GPS 
% epochs except with different reference epoch. 
% INPUT:
%   jplTime = Time used in JPL real time ephemeris [Nx1]
%
% OUTPUT:
%   epochs = GPS epoch [Nx1]
%
% See also:

epochs = jplTime+20*365*86400+0.5*86400;

end