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
        
        % Vehicle constraints
        VehicleConstraint (5)
        
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
        
        
        
    end
end

