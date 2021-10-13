classdef NavsuTests <  matlab.unittest.TestCase

    methods (Test)

        function curlDirContentTest(testCase)
            % Test the function retrieving the list of files in a cddis
            % directory.

            if ~ispc
                % check if curl exists
                [a, ~] = system('curl');
                if a  == 2
                    fileName = 'ftp://gdc.cddis.eosdis.nasa.gov/gnss/data/campaign/mgex/daily/rinex3/2020/brdm/';
                    fileList = navsu.ftp.curlGetDirectoryContents(fileName);
                    % should have exactly 176 files
                    testCase.verifyEqual(numel(fileList), 176);
                    % test the name of the last
                    testCase.verifyTrue(strcmp(fileList{end}, 'brdm1900.20p.Z'));
                else
                    warning('No curl distribution found. Cannot run this test.');
                end
            else
                warning('Could not test cURL function. Missing netrc file.');
            end

        end

        function orbitParserTest(testCase)
            % Test navigation message parsing for different sources.
            %% Parse each file
            % RINEX 2 navigation data
            fileRinex2 = fullfile(fileparts(mfilename('fullpath')), 'test-data', 'thti0500.19n');
            eph2 = navsu.readfiles.loadRinexNav(fileRinex2);

            % RINEX navigation data
            fileRinex3 = fullfile(fileparts(mfilename('fullpath')), 'test-data', 'brdm0500.19p');
            eph3 = navsu.readfiles.loadRinexNav(fileRinex3);

            % IGS final sp3
            fileIgs = fullfile(fileparts(mfilename('fullpath')), 'test-data', 'igs20412.sp3');
            ephP = navsu.readfiles.readSp3(fileIgs);

            %% Prop settings
            epochs = navsu.time.gps2epochs(2041, 180000) + (30:30:900)';
            prns = 5 * ones(size(epochs));
            constInds = 1 * ones(size(epochs));

            % Interpolate the sp3 file
            posp = navsu.geo.pephInterp(ephP, prns, constInds, epochs);

            % Propagate the broadcast navigation data from the RINEX 2
            dataB2 = navsu.geo.propNavMsg(eph2, prns, constInds, epochs);
            posb2 = [dataB2.x dataB2.y dataB2.z];

            % Propagate the broadcast navigation data from the RINEX 3
            dataB3 = navsu.geo.propNavMsg(eph3, prns, constInds, epochs);
            posb3 = [dataB3.x dataB3.y dataB3.z];

            %% Check the orbits
            % RINEX 2 and RINEX 3 data should be identical
            testCase.verifyEqual(posb2, posb3);

            % Broadcast and precise should be close
            testCase.assertLessThan(abs(posp - posb3), 1);

        end

    end

end
