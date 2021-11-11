% Unit test script to test the functions navsu.svprn.mapSignalFreq,
% navsu.svprn.prn2sv, navsu.svprn.prn2x.
% Test 2 also tests navsu.svprn.prn2FreqChanGlonass

% some shared variables here
consts = navsu.thirdparty.initConstellation(1, 1, 1, 1, 0);
% julian date for the tests corresponding to 2018/03/01
jd = 2.458179318055556e+06;



%% Test 1: GPS frequencies
cString = 'GPS';

L1 = navsu.svprn.mapSignalFreq(ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
assert(all(L1 == 1.57542e+09), ['Got ', cString, ' L1 wrong!']);

L2 = navsu.svprn.mapSignalFreq(2*ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
assert(all(L2 == 1.2276e+09), ['Got ', cString, ' L2 wrong!']);

L5 = navsu.svprn.mapSignalFreq(5*ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
assert(all(L5 == 1.17645e+09), ['Got ', cString, ' L5 wrong!']);

% now test all three at once to validate the option of running multiple
% frequencies at once

allFreq = navsu.svprn.mapSignalFreq([1 2 5] .*ones(consts.(cString).numSat, 1), ...
                                    consts.(cString).PRN, ...
                                    consts.constInds(consts.(cString).indexes));
assert(all(allFreq == [1.57542, 1.2276, 1.17645]*1e+09, 'all'), ...
       ['Got ', cString, ' triple freq. case wrong!']);

%% Test 2: GLONASS frequencies
cString = 'GLONASS';


desiredAnswer = 1.0e+09 * [1.6025625; 1.59975; 1.6048125; 1.605375; ...
                           1.6025625; 1.59975; 1.6048125; 1.605375; ...
                           1.600875; 1.5980625; 1.602; 1.6014375; ...
                           1.600875; 1.5980625; 1.602; 1.6014375; ...
                           1.60425; 1.6003125; 1.6036875; 1.603125; ...
                           1.60425; 1.6003125; 1.6036875; 1.603125];

gloL1 = navsu.svprn.mapSignalFreq(ones(consts.(cString).numSat, 1), ...
                                  consts.(cString).PRN, ...
                                  consts.constInds(consts.(cString).indexes), jd);
assert(all(gloL1 == desiredAnswer), ['Got ', cString, ' L1 wrong!']);

%% Test 3: GALILEO frequencies
cString = 'Galileo';

L1 = navsu.svprn.mapSignalFreq(ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
assert(all(L1 == 1.57542e+09), ['Got ', cString, ' E1 wrong!']);

L5 = navsu.svprn.mapSignalFreq(8*ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
assert(all(L5 == 1.191795e+09), ['Got ', cString, ' E5 wrong!']);

%% Test 4: GPS SVN numbers
cString = 'GPS';

gpsSvns = [63 61 69 NaN 50 67 48 72 68 73 46 58 43 41 55 56 53 NaN 59 ...
           51 45 47 60 65 62 71 66 44 57 64 52 70]';
activeSats = isfinite(gpsSvns);
% run with 
svn = navsu.svprn.prn2svn(consts.(cString).PRN', ...
                          jd, ...
                          consts.constInds(consts.(cString).indexes)');

% make sure we got the active ones right
assert(all(svn(activeSats) == gpsSvns(activeSats)), ...
       ['Failed to get ', cString, ' SVN numbers.']);

assert(all(isfinite(svn) == activeSats), ...
       ['Incorrect number of ', cString, ' PRNs active!']);

% now try with different input dimensions
svn = navsu.svprn.prn2svn(consts.(cString).PRN, ...
                          jd, ...
                          consts.constInds(consts.(cString).indexes)');
assert(all(svn(activeSats) == gpsSvns(activeSats)), ...
       'Failed prn2svn with column input.');

assert(all(isfinite(svn) == activeSats), ...
       'Failed prn2svn with column input.');




