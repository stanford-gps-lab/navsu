function [R, sunpos] = svLocalFrame(svPos,epochs,sunpos)
% function to create the satellite local frame rotation matrix

if nargin < 3   
    % get sun position for each time
    jd = navsu.time.epochs2jd(epochs);
    sunpos = navsu.geo.sunVecEcef(jd)';
end

% produce rotation matrices
nEpochs = length(epochs);
R = zeros(3,3,nEpochs);
for tdx = 1:nEpochs
    sunposi = sunpos(:,tdx);
    svPosi  = svPos(tdx,:)';
    
    % Build body axis rotation matrix
    e = (sunposi-svPosi)./norm(sunposi-svPosi);
    k = -svPosi./norm(svPosi);

    j = cross(k,e)/norm(cross(k,e));
    i = cross(j,k)/norm(cross(j,k));
    
    R(:,:,tdx) = [i j k];
end


end