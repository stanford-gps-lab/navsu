function [R, sunpos] = svLocalFrame(svPos, epochs, sunpos)
% function to create the satellite local frame rotation matrix
% 
% Inputs:
%  svPos    N x 3 matrix of satellite positions
% And one of the two:
%  epochs   N x 1 vector of epochs
%  sunpos   3 x N matrix of sun position at epochs

if nargin < 3 || isempty(sunpos)
    % compute sun position for each time
    jd = navsu.time.epochs2jd(epochs);
    sunpos = navsu.geo.sunVecEcef(jd)';
end

% want 3 by N sv positions
svPos = svPos';

% Build body axis rotation matrix
e = (sunpos-svPos) ./ vecnorm(sunpos-svPos, 2, 1);
k = -svPos ./ vecnorm(svPos, 2, 1);

cke = cross(k, e);
j = cke ./ vecnorm(cke, 2, 1);
cjk = cross(j, k);
i = cjk ./ vecnorm(cjk, 2, 1);

% produce rotation matrices
R = permute(cat(3, i, j, k), [1 3 2]);

end