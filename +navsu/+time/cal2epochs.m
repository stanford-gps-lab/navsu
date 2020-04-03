function epochs = cal2epochs(yr,mn,dy,hr,min,sec)

% Allow for single vector input
if nargin == 1 && size(yr,2) == 6
   mn  = yr(:,2);
   dy  = yr(:,3);
   hr  = yr(:,4);
   min = yr(:,5);
   sec = yr(:,6);
   yr  = yr(:,1);
end

jd = navsu.time.cal2jd(yr,mn,dy);
[gpsweek,sow] = navsu.time.jd2gps(jd);
epochs = gpsweek*86400*7+sow+hr*3600+min*60+sec;

end