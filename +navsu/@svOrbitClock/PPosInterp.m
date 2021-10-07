function [posP, velP, pPosInds, pPosPoly, ROut] = PPosInterp(obj, prns, ...
    constInds, epochs, pPosInds, pPosPoly, FLAG_APC_OFFSET, sunPos, dttx)

%% PPosInterp
% Given IGS precise orbit from IGS products, this function performs either
% lagrange or polynomial interpolation of the orbit to the desired epochs
%
% Required Inputs:
%  PRNs               - N-length vector of PRN per desired interpolation
%  constInds          - N-length vector of constellation index per desired
%                       interpolation.
%  epochs             - N-length GPS epoch per desired interpolation
%
% Optional Inputs:
%  pPosInds           - currently not used
%  pPosPoly           - currently not used
%  FLAG_APC_OFFSET    - flag indicating whether or not to offset the
%                       output position by the antenna phase center using 
%                       a nominal attitude model
%  atxData            - stucture containing data from IGS ATX file- can be
%                       produced by navsu.readfiles.readAtx.m
%  sunPos             - position of the sun in ECEF
%  dttx               - additional fine time offset
%
% Outputs:
%  posP               - interpolated satellite precise position
%  velP               - interpolated satellite precise velocity
%  pPosInds           - currently not used
%  pPosPoly           - currently not used
%  ROut               - rotation matrix from ECEF to SV body frame (or vice
%                       versa)
% 
% Works mostly as a wrapper for navsu.geo.pephInterp

%% check optional inputs
if nargin < 5 || isempty(pPosInds)
    pPosInds = zeros(length(epochs), 2);
else
    warning('Specified propagation indices currently not used in PposInterp')
end

if nargin < 6 || isempty(pPosPoly)
    pPosPoly = NaN(length(epochs), 6, obj.settings.polyfit.pfit+1);
end

if nargin < 7
    FLAG_APC_OFFSET = false;
end

if nargin < 8
    sunPos = [];
end

if nargin < 9 || isempty(dttx)
    dttx = zeros(length(epochs), 1);
end


% call the precise ephemeris interpolation function
[posP, velP, pPosInds, pPosPoly, ROut] = navsu.geo.pephInterp( ...
    obj.PEph, prns, constInds, epochs, ...
    'orbitInterpMethod', obj.settings.orbitInterpMethod, ...
    'nPolyFit', obj.settings.polyfit.nPolyFit, ...
    'pfit', obj.settings.polyfit.pfit, ...
    'FLAG_APC_OFFSET', FLAG_APC_OFFSET, ...
    'pPosInds', pPosInds, ...
    'pPosPoly', pPosPoly, ...
    'atxData', obj.atx, ...
    'sunPos', sunPos, ...
    'dttx', dttx);

end