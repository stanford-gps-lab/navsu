classdef MeasEnum < uint8
    % MeasEnum     enumeration of the types of measurements available
    
    enumeration
        % GNSS observation 
        GNSS (1)
        
        % Direct Position Measurement  the ECEF position in [m]
        Position (2)
        
        % Direct Velocity Measurement in ECEF in m/s
        Velocity (3)
        
        % Inertial measurements :)
        IMU (4)
        
    end
end

