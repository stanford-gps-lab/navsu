classdef DFMCnavEngineTest < matlab.unittest.TestCase
    % Unittest class for the DFMCnavEngine
    %   Runs a bank of unit tests for the DFMCnavEngine class.

    properties
        navEngine   % the nav engine object
        eph         % the ephemeris struct
    end

    methods (TestClassSetup)

        function initializeNavEngine(testCase)
            % get file path to brdc file            
            filename = fullfile(fileparts(mfilename('fullpath')), ...
                                'test-data', 'brdm0500.19p');

            testCase.eph = navsu.readfiles.loadRinexNav(filename);

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

            % make sure it worked at least for some
            testCase.verifyTrue(any(isfinite(testCase.navEngine.satPos), 'all'));
        end

    end
end