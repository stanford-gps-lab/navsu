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
    % what year should include 
    ephRDate = Eph_R(2:7,1)';
    ephRDate(1) = ephRDate(1)+2000;
    
    [~,~,~,leapSecond] = navsu.time.utc2gps(ephRDate,1);
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
end



















