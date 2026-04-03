% Unit test script to test the functions navsu.svprn.mapSignalFreq,
% navsu.svprn.prn2sv, navsu.svprn.prn2x.
% Test 2 also tests navsu.svprn.prn2FreqChanGlonass

classdef (SharedTestFixtures={ matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath'))))}) ...
        svprnTest < matlab.unittest.TestCase

properties
consts = navsu.thirdparty.initConstellation(1, 1, 1, 1, 0);
% julian date for the tests corresponding to 2018/03/01
jd = 2.458179318055556e+06;
end

methods (Test)

function prn2svTest(tc)


%% Test 1: GPS frequencies
cString = 'GPS';
consts = tc.consts;

L1 = navsu.svprn.mapSignalFreq(ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
tc.verifyEqual(size(L1), [32, 1])
tc.verifyEqual(L1, 1.57542e+09 * ones(size(L1)), ['Got ', cString, ' L1 wrong!'])

L2 = navsu.svprn.mapSignalFreq(2*ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
tc.verifyEqual(size(L2), [32, 1])
tc.verifyEqual(L2, 1.2276e+09* ones(size(L2)), ['Got ', cString, ' L2 wrong!']);

L5 = navsu.svprn.mapSignalFreq(5*ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
tc.verifyEqual(size(L5), [32, 1])
tc.verifyEqual(L5, 1.17645e+09* ones(size(L5)) , ['Got ', cString, ' L5 wrong!']);

% now test all three at once to validate the option of running multiple
% frequencies at once

allFreq = navsu.svprn.mapSignalFreq([1 2 5] .*ones(consts.(cString).numSat, 1), ...
                                    consts.(cString).PRN, ...
                                    consts.constInds(consts.(cString).indexes));
tc.verifyEqual(allFreq, repmat([1.57542, 1.2276, 1.17645]*1e+09, [32,1]), ...
       ['Got ', cString, ' triple freq. case wrong!'])
end

function TestGlonass(tc)
%% Test 2: GLONASS frequencies
cString = 'GLONASS';
consts = tc.consts;

% just to G1
desiredAnswerG1 = 1.0e+09 * [1.6025625; 1.59975; 1.6048125; 1.605375; ...
                           1.6025625; 1.59975; 1.6048125; 1.605375; ...
                           1.600875; 1.5980625; 1.602; 1.6014375; ...
                           1.600875; 1.5980625; 1.602; 1.6014375; ...
                           1.60425; 1.6003125; 1.6036875; 1.603125; ...
                           1.60425; 1.6003125; 1.6036875; 1.603125];

gloG1 = navsu.svprn.mapSignalFreq(ones(consts.(cString).numSat, 1), ...
                                  consts.(cString).PRN, ...
                                  consts.constInds(consts.(cString).indexes), tc.jd);
tc.verifyEqual(gloG1, desiredAnswerG1, ['Got ', cString, ' G1 wrong!']);

% try G1, G2, G2a at once plus an illegal one
desiredAnswerG12NaN3 = 1e6*[1602, 1246, NaN, 1248.06] ...
              + 1e6*[9/16, 7/16, 0, 0].*navsu.svprn.prn2FreqChanGlonass(consts.(cString).PRN, tc.jd);
gloAll3 = navsu.svprn.mapSignalFreq([1 2 5 6] .* ones(consts.(cString).numSat, 1), ...
                                    consts.(cString).PRN, ...
                                    consts.constInds(consts.(cString).indexes), tc.jd);
freqExists = isfinite(desiredAnswerG12NaN3);
tc.verifyEqual(gloAll3(freqExists), desiredAnswerG12NaN3(freqExists), ...
    ['Got ', cString, ' triple freq. wrong!']);
tc.verifyEqual(freqExists, isfinite(gloAll3), ...
    ['Got ', cString, ' number of legal GLONASS signals wrong!']);

end

function TestGalileo(tc)
%% Test 3: GALILEO frequencies
cString = 'Galileo';
consts = tc.consts;
E1 = navsu.svprn.mapSignalFreq(ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
tc.verifyEqual(E1, 1.57542e+09 * ones(size(E1)), ['Got ', cString, ' E1 wrong!']);

E5 = navsu.svprn.mapSignalFreq(8*ones(consts.(cString).numSat, 1), ...
                               consts.(cString).PRN, ...
                               consts.constInds(consts.(cString).indexes));
tc.verifyEqual(E5, 1.191795e+09* ones(size(E1)), ['Got ', cString, ' E5 wrong!']);
end

function TestMultiFreq(tc)
%% Test 4: Multi-frequency frequencies
consts = tc.consts;

L1G2EillegalE5a = navsu.svprn.mapSignalFreq( ...
    [[1 0] .* ones(consts.GPS.numSat, 1); ...
     [2 0] .* ones(consts.GLONASS.numSat, 1); ...
     [3 5] .* ones(consts.Galileo.numSat, 1)], ...
    [consts.GPS.PRN'; consts.GLONASS.PRN'; consts.Galileo.PRN'], ...
    consts.constInds([consts.GPS.indexes, consts.GLONASS.indexes, consts.Galileo.indexes])', ...
    tc.jd);

expectedOutput = [[1575.42e+06 NaN].*ones(consts.GPS.numSat, 1); ...
                  [1246e+06 NaN]+7/16*1e+06*navsu.svprn.prn2FreqChanGlonass(consts.GLONASS.PRN, tc.jd); ...
                  [NaN 1176.45e+06].*ones(consts.Galileo.numSat, 1)];

freqExists = isfinite(expectedOutput);
tc.verifyEqual(L1G2EillegalE5a(freqExists), expectedOutput(freqExists), ...
       'Got multi constellation frequencies wrong!');
tc.verifyEqual(freqExists, isfinite(L1G2EillegalE5a), ...
       'Got number of legal multi-frequency signals wrong!');
end

%% Test 5: GPS SVN numbers
function TestSVN(tc)
cString = 'GPS';
consts = tc.consts;

gpsSvns = [63 61 69 NaN 50 67 48 72 68 73 46 58 43 41 55 56 53 NaN 59 ...
           51 45 47 60 65 62 71 66 44 57 64 52 70]';
activeSats = isfinite(gpsSvns);

svn = navsu.svprn.prn2svn(consts.(cString).PRN', ...
                          tc.jd, ...
                          consts.constInds(consts.(cString).indexes)');

% make sure we got the active ones right
tc.verifyEqual(svn(activeSats), gpsSvns(activeSats), ...
       ['Failed to get ', cString, ' SVN numbers.']);

tc.verifyEqual(isfinite(svn), activeSats, ...
       ['Incorrect number of ', cString, ' PRNs active!']);

% now try with different input dimensions
svn = navsu.svprn.prn2svn(consts.(cString).PRN, ...
                          tc.jd, ...
                          consts.constInds(consts.(cString).indexes)');
tc.verifyEqual(svn(activeSats), gpsSvns(activeSats), ...
       'Failed prn2svn with column input.');

tc.verifyEqual(isfinite(svn), activeSats, ...
       'Failed prn2svn with column input.');
end

function TestConstellationDefault(tc)
consts = tc.consts;
%% Test 6: Default to constellation GPS
svn = navsu.svprn.prn2svn(consts.GPS.PRN, tc.jd);
gpsSvns = [63 61 69 NaN 50 67 48 72 68 73 46 58 43 41 55 56 53 NaN 59 ...
           51 45 47 60 65 62 71 66 44 57 64 52 70]';
activeSats = isfinite(gpsSvns);

tc.verifyEqual(svn(activeSats), gpsSvns(activeSats), ...
       'prn2svn failed to default to GPS.');

tc.verifyEqual(isfinite(svn), activeSats, ...
       'prn2svn failed to default to GPS.');
end

function TestPRN2SVN(tc)
consts = tc.consts;
%% Test 7: PRN 2 SVN 2 PRN
useConsts = {'GPS'; 'GLONASS'; 'Galileo'; 'BeiDou'};
prnCell = cellfun(@(x) consts.(x).PRN, useConsts, 'UniformOutput', false);
PRNs = horzcat(prnCell{:})';
constCell = cellfun(@(x) consts.(x).indexes, useConsts, 'UniformOutput', false);
constIds = consts.constInds(horzcat(constCell{:})');
prnAfter = navsu.svprn.svn2prn(navsu.svprn.prn2svn(PRNs, tc.jd, constIds), ...
                               tc.jd, ...
                               constIds);
prnExists = isfinite(prnAfter);
tc.verifyEqual(prnAfter(prnExists), PRNs(prnExists), ...
       'Converting PRNs to SVNs and back failed.');
end

end

end