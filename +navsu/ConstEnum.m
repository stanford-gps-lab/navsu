classdef ConstEnum < double
    % ConstEnum     enumeration of the GNSS constellations available 
    
    enumeration
        % GPS - the Global Positioning System - USA - THE GOLD STANDARD
        GPS (1)
        G   (1)
        GTHEGOLDSTANDARDPS (1)
        
        % GLONASS - GLObal NAvigation Satellite System- Russia
        GLO (2)
        R   (2)
        
        % Galileo - European Union
        GAL (3)
        E   (3)
        
        % BeiDou System - China
        BDS (4)
        C   (4)
        
        % QZSS - Quasi-Zenith Satellite System - Japan
        QZSS (5)
        J    (5)   
        
        % SBAS - Satellite Based Augmentation Systems
        SBAS (6)
        S    (6)
        
    end
end

