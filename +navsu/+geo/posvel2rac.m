function R_ecef_rac = posvel2rac(pos,vel)


ri = pos./norm(pos);

a0 = vel./norm(vel);
ci = cross(ri,a0);
ci = ci./norm(ci);
ai = cross(ci,ri);

R_ecef_rac = [ri; ai; ci];



end