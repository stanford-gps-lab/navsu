function Beph = multiGNSSEphArray2Struct(BephA,filename)
% 
% Beph.gps = [];
% Beph.glo = [];
% Beph.gal = [];
% Beph.bds = [];
% Beph.qzss = [];
% Beph.iono = [];
% Beph.leapSecond = [];

Beph.gps  = navsu.readfiles.ephArray2Struct(BephA.gps',filename,BephA.leapSecond,'GPS');
Beph.glo  = navsu.readfiles.ephArray2StructGlonass(BephA.glo',filename,BephA.leapSecond);
Beph.gal  = navsu.readfiles.ephArray2Struct(BephA.gal',filename,BephA.leapSecond,'GAL');
Beph.bds  = navsu.readfiles.ephArray2Struct(BephA.bds',filename,BephA.leapSecond,'BDS');
Beph.qzss = navsu.readfiles.ephArray2Struct(BephA.qzss',filename,BephA.leapSecond,'QZSS');
Beph.iono = BephA.iono;

end