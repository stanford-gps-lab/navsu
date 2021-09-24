function [phwindup] = carrierPhaseWindupGGM(epochi, XR, XS, phwindup)

% SYNTAX:
%   [phwindup] = phase_windup_correction(time, XR, XS, SP3, phwindup);
%
% INPUT:
%   time = GPS time
%   XR   = XYZ receiver position: N by 3 matrix or 1 by 3 vector
%   XS   = XYZ satellite position N by 3 matrix
%   phwindup = phase wind-up (previous value) N by 1 vector
%
% OUTPUT:
%   phwindup = phase wind-up (updated value)
%
% DESCRIPTION:
%   Computation of the phase wind-up terms.

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 3
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
%   
%   Completely rewritten for vectorization by Fabian Rothmaier, 09/2021

% redo vectorized
llh = navsu.geo.xyz2llh(XR);
phi = llh(:, 1)'*pi/180;
lam = llh(:, 2)'*pi/180;

% east (a) and north (b) local unit vectors
a = [-sin(lam); cos(lam); zeros(size(lam))];
b = [-sin(phi).*cos(lam); -sin(phi).*sin(lam); cos(phi)];

% satellite-fixed local unit vectors
R = navsu.geo.svLocalFrame(XS, epochi);
i = squeeze(R(:, 1, :));
j = squeeze(R(:, 2, :));
k = squeeze(R(:, 3, :));

%receiver and satellites effective dipole vectors
nSat = size(XS, 1);
Dr = a - k.*dot(k, a.*ones(1, nSat)) + cross(k, b.*ones(1, nSat));
Ds = i - k.*dot(k, i) - cross(k, j);

%phase wind-up computation
psi = dot(k, cross(Ds, Dr));
arg = dot(Ds, Dr) ./ (vecnorm(Ds, 2, 1) .* vecnorm(Dr, 2, 1));
% limit to [-1 1]
arg = min(max(arg, -1), 1);
dPhi = (sign(psi) .* acos(arg) / (2*pi))';

% add up with previous values
N = round(phwindup - dPhi);
N(phwindup == 0) = 0;

phwindup = dPhi + N;
phwindup(isnan(phwindup)) = 0;

end
