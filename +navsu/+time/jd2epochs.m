function epochs = jd2epochs(jds)


[weeks,tow] = navsu.time.jd2gps(jds);

epochs = weeks*604800+tow;



end