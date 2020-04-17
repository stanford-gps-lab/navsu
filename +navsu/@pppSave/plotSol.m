function plotSol(obj,varargin)

% 
p = inputParser;

p.addParameter('truePosEcef',[]);

% parse the results
parse(p, varargin{:});
res = p.Results;
truePosEcef          = res.truePosEcef;          % truth position to compare to

% Plot the position and clock bias in ENU

xyz  = obj.posSave;
b    = obj.clockSave;
epochs = obj.epochSave;
llh0 = navsu.geo.xyz2llh(xyz(1,:));

utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));

enu = nan(size(xyz));
for idx = 1:size(xyz,1)
    [enu(idx,1),enu(idx,2),enu(idx,3)] = ...
        navsu.thirdparty.cart2utm(xyz(idx,1),xyz(idx,2),xyz(idx,3),utmZone);
end

%
if ~isempty(truePosEcef) 
    truePosEnu = nan(1,3);
    [truePosEnu(1), truePosEnu(2), truePosEnu(3)] = navsu.thirdparty.cart2utm(truePosEcef(1),truePosEcef(2),truePosEcef(3),utmZone);
end


%% plot
figure; 
ha = navsu.thirdparty.tightSubplot(4,1,0.02,[0.1 0.1],[0.07 0.05]);

if ~isempty(truePosEnu)
    compPosEnu = truePosEnu;
else
   compPosEnu = enu(end,:); 
end

yplot = [enu-compPosEnu b*navsu.constants.c]';
tplot = (epochs-epochs(1))/60; % minutes
ylabels = {'East [m]' 'North [m]' 'Up [m]' 'Clock [m]'};
for idx = 1:4
    axes(ha(idx))
    plot(tplot,yplot(idx,:))
    ylabel(ylabels{idx})
    grid on
    if idx < 4
       xticklabels('') 
    else
       xlabel('Time [min]') 
    end
end

end