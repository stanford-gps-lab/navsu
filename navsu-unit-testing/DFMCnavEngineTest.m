classdef (SharedTestFixtures={ matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename('fullpath'))))}) ...
        DFMCnavEngineTest < matlab.unittest.TestCase
    % Unittest class for the DFMCnavEngine
    %   Runs a bank of unit tests for the DFMCnavEngine class.

    properties
        navEngine   % the nav engine object
        eph         % the broadcast ephemeris struct
    end

    methods (TestClassSetup)

        function initializeNavEngine(testCase)
            % get file path to brdc file
            brdcFilename = fullfile(fileparts(mfilename('fullpath')), ...
                                    'test-data', 'brdm0500.19p');

            testCase.eph = navsu.readfiles.loadRinexNav(brdcFilename);

            testCase.navEngine = navsu.lsNav.DFMCnavigationEngine(testCase.eph);

        end

    end

    % this would be needed if e.g. a new brdc file was downloaded
    %     methods(TestClassTeardown)
    %
    %         function deleteDownloadedFiles(testCase)
    %
    %         end
    %
    %     end

    methods (Test)

        function testOrbitProp(testCase)
            % Test the orbit propagation

            epProp = navsu.time.gps2epochs(testCase.eph.gps.GPS_week_num, ...
                                           testCase.eph.gps.Toe);

            % attempt orbit propagation for all satellites
            testCase.navEngine.propagateOrbits(1:testCase.navEngine.numSats, ...
                                               mean(epProp));

            % make sure it we got results at least for some
            testCase.verifyTrue(any(isfinite(testCase.navEngine.satPos), 'all'));
        end
        
        function testRinexParser(testCase)
            
            % get path to Rinex file
            navsuPath = fileparts(fileparts(mfilename('fullpath')));
            obsFile = fullfile(navsuPath, 'examples', 'swift-gnss-20200312-093212.obs');
            testCase.assumeThat(obsFile, matlab.unittest.constraints.IsFile);
            % read obs file
            [obsStruc, constellations, epochs] = navsu.readfiles.loadRinexObs(obsFile);
            obsGnssRaw.meas      = obsStruc;
            obsGnssRaw.PRN       = constellations.PRN;
            obsGnssRaw.constInds = constellations.constInds;
            obsGnssRaw.epochs    = epochs;
            
            % parse data for some epoch for all satellites
            ep = 100;
            obsData = testCase.navEngine.readRinexData(obsGnssRaw, [], ep);
            
            % make sure every satellite is represented
            testCase.verifyEqual(size(obsData.code, 1), ...
                                 size(obsGnssRaw.meas.C1C, 1));
                             
            % make sure the code values match
            fn = fieldnames(obsGnssRaw.meas);
            codeData = cellfun(@(x) strcmp(x(1), 'C') && length(x) == 3, fn);
            
            for fni = find(codeData(:))'
                % check every code measurement
                rnxData = obsGnssRaw.meas.(fn{fni})(:, ep);
                haveMeas = isfinite(rnxData) & rnxData ~= 0;
                
                testCase.verifyTrue(all(any(obsData.code(haveMeas, :) ...
                                            == rnxData(haveMeas), 2)))
            end
            
        end

    end
end
