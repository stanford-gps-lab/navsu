function cbias = clockInterp(obj, prns, constInds, epochs)
% clockInterp
% DESCRIPTION:
% Very simple linear interpolation of IGS clock products.  This is called
% by navsu.svOrbitClock.clock
%
% INPUTS:
%  prns                - N-length vector of PRN per desired interpolation
%  constInds           - N-length vector of constellation index per desired
%                        interpolation
%  epochs              - N-length GPS epoch per desired interpolation
%  Clck                - IGS precise clock structure output from loadCFst.m
%
% OUTPUTS:
%  cbias               - N-length vector of interpolated clock bias per
%                        desired interpolation
%
% See also: navsu.svOrbitClock.clock
%%

% retrieve used object property
Clck = obj.PClock;

cbias = nan(size(prns));
cbias = cbias(:);
for idx = 1:length(prns)
   satInd = find(Clck.PRNs == prns(idx) & Clck.constInds == constInds(idx));   
   eInd1 = find(Clck.Cepochs <= epochs(idx), 1, 'last'); 
   eInd2 = find(Clck.Cepochs >  epochs(idx), 1, 'first'); 

   if isempty(satInd) || isempty(eInd1) || isempty(eInd2)
       continue;
   end

   % Simple linear interpolation between points before and after
   cbias(idx) = (Clck.Cclk(satInd,eInd2)-Clck.Cclk(satInd,eInd1)) * (epochs(idx)-Clck.Cepochs(eInd1)) ./ ...
                (Clck.Cepochs(eInd2)-Clck.Cepochs(eInd1)) + Clck.Cclk(satInd,eInd1);
%    cbias(idx) = interp1(Clck.Cepochs([eInd1 eInd2]),Clck.Cclk(satInd,[eInd1 eInd2]),epochs(idx));
   
end


end