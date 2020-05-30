classdef MeasEnum < uint8
    % MeasEnum     enumeration of the types of measurements available
    
    enumeration
        % Any 
        AnyMeas (0)
        
        % GNSS observation 
        GNSS (1)
        
        % Direct Position Measurement  the ECEF position in [m]
        Position (2)
        
        % Direct Velocity Measurement in ECEF in m/s
        Velocity (3)
        
        % Inertial measurements :)
        IMU (4)
        
        % Vehicle constraints
        VehicleConstraint (5)
        
        % Vehicle data measurement - wheel speed, steering angle 
        Wheels (6)
        
        %% Measurement subtypes- 
        % GNSS Code phase
        Code (10)
        
        % GNSS Carrier phase
        Carrier (11)
        
        % GNSS Doppler
        Doppler (12)
        
        %% Position subtypes
        % ECEF X position
        PosX (21)
        
        % ECEF Y position
        PosY (22)
        
        % ECEF Z position
        PosZ (23)
        
        %% Velocity subtype
        % ECEF X position
        VelX (31)
        
        % ECEF Y position
        VelY (32)
        
        % ECEF Z position
        VelZ (33)
        
        %% Vehicle contraint subtypes
        % Vertical constraint
        NoSlipVertical (51)
        
        % Side slip contraint
        NoSlipCross (52)
        
        %% Vehicle speed measurements
        % Front left wheel speed 
        SpeedFrontLeft (60)
        
        % Front right wheel speed
        SpeedFrontRight (61)
        
        % Rear left wheel speed
        SpeedRearLeft (62)
        
        % Rear right wheel speed
        SpeedRearRight (63)
        
        % Steering angle
        SteeringAngle (64)
        
        % Vehicle speed
        VehicleSpeed (65)
        
        % Stationary flag
        Stationary (66)
        
        % Transmission setting
        TranmissionSetting (67)
        
        % Parking brake status
        ParkingBrakeStatus (68)
        
    end
end

