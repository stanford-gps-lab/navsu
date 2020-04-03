function epochs = gps2epochs(gpsWeek,tow)

% gpsWeek = floor(epochs/(604800));
% tow = mod(epochs,604800);

epochs = gpsWeek*604800+tow;

end