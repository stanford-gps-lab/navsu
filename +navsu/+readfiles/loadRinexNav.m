function eph = loadRinexNav(filename, varargin)
% loadRinexNav
% DESCRIPTION:
%   Parse a multi (or single) constellation RINEX navigation message file
% INPUT:
%   filename       = name of the RINEX navigation file with the path
%   
% OPTIONAL INPUTS:
%   outFormat      = 'struct' or 'array' - specifies the format of the output.
%                    Default is 'struct'
%   constellations = constellation object created by
%                    navsu.readfiles.initConstellation and can be used to
%                    select what constellations are parsed. Default parses
%                    everything available. 
%
% OUTPUT:
%   eph     = structure containing navigation message data from each
%             constellation.  Default output forms each constellation
%             output as a structure, but can be assigned to be arrays using
%             the 'outFormat' input
%             Units of each element should be the same as what was in the
%             RINEX file. 
%
% See also: navsu.readfiles.loadRinexObs
%
% This has been heavily modified from the original goGPS code.
%
%--------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2013 Mirko Reguzzoni,Eugenio Realini
% Portions of code contributed by Damiano Triglione (2012)
%--------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

p = inputParser;
p.addParameter('outFormat', 'struct'); % 'struct' or 'array'
p.addParameter('constellations', navsu.readfiles.initConstellation(1,1,1,1,1,1));

% parse the results
parse(p, varargin{:});
res = p.Results;
constellations = res.constellations;
outFormat      = res.outFormat;

Eph_G = []; 
Eph_R = [];
Eph_E = [];
Eph_C = [];
Eph_J = [];
iono = zeros(8,1);

if constellations.GPS.enabled
    %parse RINEX navigation file (GPS) NOTE: filename expected to
    %end with 'n' or 'N' (GPS) or with 'p' or 'P' (mixed GNSS)
    constellationsi = navsu.readfiles.initConstellation(1,0,0,0,0,0);
    
    [Eph_G, iono,glut,leapSecond] = navsu.readfiles.rinexGetNav(filename, constellationsi);
end

if constellations.GLONASS.enabled
    %parse RINEX navigation file (GLONASS)
    constellationsi = navsu.readfiles.initConstellation(0,1,0,0,0,0);
    
    [Eph_R, iono,glut,leapSecond] = navsu.readfiles.rinexGetNav(filename, constellationsi);
end

if constellations.Galileo.enabled
    %parse RINEX navigation file (Galileo)
    constellationsi = navsu.readfiles.initConstellation(0,0,1,0,0,0);
    
    [Eph_E, iono,glut,leapSecond] = navsu.readfiles.rinexGetNav(filename, constellationsi);
end

if constellations.BeiDou.enabled
    constellationsi = navsu.readfiles.initConstellation(0,0,0,1,0,0);
    
    %parse RINEX navigation file (BeiDou)
    [Eph_C, iono,glut,leapSecond] = navsu.readfiles.rinexGetNav(filename, constellationsi);
end

if constellations.QZSS.enabled
    constellationsi = navsu.readfiles.initConstellation(0,0,0,0,1,0);
    
    %parse RINEX navigation file (QZSS)
    [Eph_J, iono,glut,leapSecond] = navsu.readfiles.rinexGetNav(filename, constellationsi);
end

if any(isnan(leapSecond)) && ~isempty(Eph_R)
    % Pull the actual leap second count
    dateVec = Eph_R(2:7,1)';
    % Solve year ambiguity. See RINEX 3 section 6.10
    if dateVec(1) < 80
        dateVec(1) = dateVec(1) + 2000;
    elseif dateVec(1) <= 99
        dateVec(1) = dateVec(1) + 1900;
    end
    [~,~,~,leapSecond] = navsu.time.utc2gps(dateVec,1);
end

% Collect everything for output
if strcmp(outFormat,'array')
    eph.gps = Eph_G;
    eph.glo = Eph_R;
    eph.gal = Eph_E;
    eph.bds = Eph_C;
    eph.qzss = Eph_J;
    eph.iono = iono;
    eph.leapSecond = leapSecond;
elseif strcmp(outFormat,'struct')
    eph.gps  = navsu.readfiles.ephArray2Struct(Eph_G',filename,leapSecond,'GPS');
    eph.glo  = navsu.readfiles.ephArray2StructGlonass(Eph_R',filename,leapSecond);
    eph.gal  = navsu.readfiles.ephArray2Struct(Eph_E',filename,leapSecond,'GAL');
    eph.bds  = navsu.readfiles.ephArray2Struct(Eph_C',filename,leapSecond,'BDS');
    eph.qzss = navsu.readfiles.ephArray2Struct(Eph_J',filename,leapSecond,'QZSS');
    eph.iono = iono;

    % add iono correction to each constellation
    fn = fieldnames(eph);

    for fni = 1:length(fn)
        if numel(eph.(fn{fni})) > 1 || ~isstruct(eph.(fn{fni}))
            % skip "iono" field which is a struct and consts for which I
            % don't have ephemeris data
            continue
        end

        % set default values for this const
        eph.(fn{fni}).ionoCorrCoeffs = NaN(1, 8);
        
        if isstruct(iono) && isfield(iono, 'ionoCorrType')
            % do we have iono correction values for this constellation?
            constMatch = arrayfun(@(x) startsWith(fn{fni}, lower(x.ionoCorrType)), iono);
    
            if sum(constMatch) > 1
                warning(['Can not uniquely identify ', fn{fni}, ' iono correction.', ...
                    ' Not adding it to specific eph constellation struct.']);
            elseif any(constMatch)
                eph.(fn{fni}).ionoCorrCoeffs = iono(constMatch).ionoCorrCoeffs';
            elseif numel(iono) == 1 && strcmp(iono.ionoCorrType, 'RINEX2_A0-B3')
                % check for Rinex 2 case
                eph.gps.ionoCorrCoeffs = iono.ionoCorrCoeffs';
            end
        end
    end

    % also add year, doy
    jds = navsu.time.gps2jd(eph.gps.GPS_week_num, eph.gps.Toe);
    [doy, eph.year] = navsu.time.jd2doy(mean(jds));
    eph.doy = floor(doy);
end



















