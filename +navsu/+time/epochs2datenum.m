function tDatenum = epochs2datenum(epochs)
% Convert from GPS epoch to MATLAB datenum

tDatenum = datenum(navsu.time.epochs2cal(epochs,1));


end