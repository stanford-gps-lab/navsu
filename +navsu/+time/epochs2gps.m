function [gpsWeek, gpsTow] = epochs2gps(epochs)

gpsWeek = floor(epochs./(7*86400));
gpsTow = epochs-gpsWeek*7*86400;


end