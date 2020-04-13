function [R, sunpos] = svLocalFrame(svPos,epochs,sunpos)
% function to create the satellite local frame rotation matrix
% offsetECEF = R(:,:,i)*offset;

if nargin < 3
%     if isempty(strfind(path,'\mice\src\mic'))
%         AttachToMice();
%     end
    
    % get sun position for each time
    jd = utility.time.epochs2jd(epochs);
    sunpos = zeros(3,length(epochs));
    for i = 1:length(jd)
        et          = cspice_str2et(['jd ' num2str(jd(i))]);
        sunposi     = cspice_spkezr( 'sun',et , 'itrf93', 'none', 'earth');
        sunpos(:,i) = sunposi(1:3)*1000;
    end
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
    % yhat = cross(e,k)./norm(cross(e,k));
    j = cross(k,e)/norm(cross(k,e));
    i = cross(j,k)/norm(cross(j,k));
    
    R(:,:,tdx) = [i j k];
end




end