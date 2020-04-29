classdef MeasEnum
    % MeasEnum     enumeration of the types of measurements available
    
    enumeration
        % GNSS observation 
        GNSS
        
        % Direct Position Measurement  the ECEF position in [m]
        Position
        
        % Direct Velocity Measurement in ECEF in m/s
        Velocity
        
        % Inertial measurements :)
        IMU
        
    end
end

