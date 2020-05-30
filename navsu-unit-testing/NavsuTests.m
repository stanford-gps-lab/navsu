classdef NavsuTests <  matlab.unittest.TestCase
    
    
    
    methods (Test)
        function rinex2NavParserTest(testCase)
            % Test RINEX 2 navigation message parsing
            fileRinex2 = fullfile(fileparts( mfilename('fullpath')), 'data', 'thti0500.19n');
            
            eph = navsu.readfiles.loadRinexNav(fileRinex2);
            
            % Propagate the data
            
        end
        
        function rinex3NavParserTest(testCase)
            % Test RINEX 3 navigation message parsing
            fileRinex3 = fullfile(fileparts( mfilename('fullpath')), 'data', 'brdm0500.19p');
            
            eph = navsu.readfiles.loadRinexNav(fileRinex3);
            
        end
    end
    
    
    
    
    
end