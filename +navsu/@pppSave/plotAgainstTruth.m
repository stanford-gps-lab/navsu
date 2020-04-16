function [epochs,posTrue,posEst,pl,errEnu,plLoc] = plotAgainstTruth(outStruc,filenameTruth,PARAMS,varargin)

%% Parse inputs
p = inputParser;
p.addParameter('plot_flag',true);
p.addParameter('tplot_mode','time'); % alternatively, 'ind'
p.addParameter('plot_all',false);
p.addParameter('plot_pl',true);
p.addParameter('plot_google_earth',false);
p.addParameter('plot_vel_enu',false);
p.addParameter('plot_imu_bias',false);
p.addParameter('plot_pos_err_body_frame',false);
p.addParameter('plot_google_maps',false);
p.addParameter('plot_speed',false);
p.addParameter('plot_ground_track',false);
p.addParameter('plot_att',false);

% parse the results
parse(p, varargin{:});
res        = p.Results;
plot_all               = res.plot_all;
tplot_mode              = res.tplot_mode;
plot_pl                 = res.plot_pl;
plot_google_earth       = res.plot_google_earth;
plot_vel_enu            = res.plot_vel_enu;
plot_imu_bias           = res.plot_imu_bias;
plot_pos_err_body_frame = res.plot_pos_err_body_frame;
plot_google_maps        = res.plot_google_maps;
plot_speed              = res.plot_speed;
plot_ground_track       = res.plot_ground_track;
plot_att                = res.plot_att;

%%
% parse the truth data
if ~exist('posTruth','var')
    [~,posTruth, epochsTruth,cov,~,~,velEnu,stdEnu] = parsePosSolFile(filenameTruth);
    
    posMeas.pos = posTruth;
    posMeas.epochs = epochsTruth;
    
    posMeas.enuPos = nan(size(posTruth));
    posMeas.enuStd = stdEnu;
    
    llh0 = xyz2llh(posTruth(1,:));
    
    utmZone = findUtmZone(llh0(1),llh0(2));
    for idx = 1:size(posTruth,1)
        [posMeas.enuPos(idx,1),posMeas.enuPos(idx,2),posMeas.enuPos(idx,3)] = ...
            cart2utm(posTruth(idx,1),posTruth(idx,2),posTruth(idx,3),utmZone);
    end
end

%% Pull data out of the output saving object
pos = outStruc.posSave;
vel = outStruc.velSave;
att = outStruc.attSave;
attEnu = outStruc.attEnuSave;
biasSave = outStruc.biasSave;
epochSave = outStruc.epochSave;
clockSave = outStruc.clockSave;
covEnu = outStruc.covEnuFull';

%% Compare to truth data
if ~isempty(posMeas)
    % compute the output position to ENU
    posOutEnu = nan(size(pos));
    for idx = 1:size(pos,1)
        [posOutEnu(idx,1),posOutEnu(idx,2),posOutEnu(idx,3)] = ...
            cart2utm(pos(idx,1),pos(idx,2),pos(idx,3),utmZone);
    end
    
    % Convert velocity to ENU velocity
    velOutEnu = nan(size(pos));
    for idx = 1:size(pos,1)
        velOutEnu(idx,:) = XYZ2ENU(vel(idx,:),llh0(1)*pi/180,llh0(2)*pi/180);
    end
    
    % Interpolate truth positions to compare
    truthEnuInterp = nan(size(posOutEnu));
    truthXyzInterp = nan(size(posOutEnu));
    if isempty(posMeas.epochs)
        truthEnuInterp = repmat(posMeas.enuPos,size(posOutEnu,1),1);
        truthXyzInterp = repmat(posMeas.pos,size(posOutEnu,1),1);
    else
        % Pad truth gaps with NaNs
        indsGap = find(diff(posMeas.epochs) > 1);
        posMeas.enuPos(indsGap,:) = nan(length(indsGap),3);
        posMeas.pos(indsGap,:) = nan(length(indsGap),3);
        
        interpMethod = 'linear';
        truthEnuInterp(:,1) = interp1(posMeas.epochs,posMeas.enuPos(:,1),epochSave,interpMethod);
        truthEnuInterp(:,2) = interp1(posMeas.epochs,posMeas.enuPos(:,2),epochSave,interpMethod);
        truthEnuInterp(:,3) = interp1(posMeas.epochs,posMeas.enuPos(:,3),epochSave,interpMethod);
        
        truthXyzInterp(:,1) = interp1(posMeas.epochs,posMeas.pos(:,1),epochSave,interpMethod);
        truthXyzInterp(:,2) = interp1(posMeas.epochs,posMeas.pos(:,2),epochSave,interpMethod);
        truthXyzInterp(:,3) = interp1(posMeas.epochs,posMeas.pos(:,3),epochSave,interpMethod);
    end
    
    %     velOutEnu = [nan nan nan; diff(posOutEnu)./diff(epochSave)];
    if ~isempty(velEnu)
        truthVelEnu = nan(size(pos));
        truthVelEnu(:,1) = interp1(posMeas.epochs,velEnu(:,1),epochSave,interpMethod);
        truthVelEnu(:,2) = interp1(posMeas.epochs,velEnu(:,2),epochSave,interpMethod);
        truthVelEnu(:,3) = interp1(posMeas.epochs,velEnu(:,3),epochSave,interpMethod);
        
    else
        truthVelEnu = [nan nan nan; diff(truthEnuInterp)./diff(epochSave)];
        truthVelEnu = [nan nan nan; (truthEnuInterp(3:end,:)-truthEnuInterp(1:(end-2),:))./(epochSave(3:end)-epochSave(1:(end-2))); ...
            nan nan nan];
    end
    
    if ~isempty(stdEnu)
        indsGap = find(diff(posMeas.epochs) > 1);
        posMeas.enuPos(indsGap,:) = nan(length(indsGap),3);
        
        interpMethod = 'linear';
        stdEnuInterp(:,1) = interp1(posMeas.epochs,stdEnu(:,1),epochSave,interpMethod);
        stdEnuInterp(:,2) = interp1(posMeas.epochs,stdEnu(:,2),epochSave,interpMethod);
        stdEnuInterp(:,3) = interp1(posMeas.epochs,stdEnu(:,3),epochSave,interpMethod);
    end
    
    truthAttEnu = nan(size(truthVelEnu));
    for idx = 1:size(truthVelEnu,1)
        truthAttEnu(idx,:) = initializeAttitude(truthVelEnu(idx,:)',PARAMS);
    end
    
    errEnu = posOutEnu-truthEnuInterp;
    errVelEnu = velOutEnu -truthVelEnu;
    %
    if strcmp(tplot_mode,'time')
        epochsPlot = (epochSave-epochSave(1))/60;
    else
        epochsPlot = 1:length(epochSave);
    end
    
    if 0
        indsGnss = find(epochSave == floor(epochSave));
    else
        indsGnss = find(epochSave == epochSave);
    end
    
    if plot_pl || plot_all
        fig1 = figure; clf;
        ha = tight_subplot(3,1,0.05,[0.1 0.1],[0.07 0.05]);
        labels = {'East [m]','North [m]','Up [m]'};
        for idx = 1:3
            axes(ha(idx))
            
            area(epochsPlot(indsGnss),sqrt(covEnu(indsGnss,idx))*3)
            hold on
            if ~isempty(stdEnu)
                % Also shade the covariance of the truth estimate
                area(epochsPlot(indsGnss),stdEnuInterp(indsGnss,idx),'FaceColor','y')
            end
            hold on;
            plot(epochsPlot(indsGnss),abs(errEnu(indsGnss,idx)),'.');
            
            
            
            
            %         semilogy(epochsPlot(indsGnss),sqrt(covEnu(indsGnss,idx))*3);
            plot(epochsPlot(indsGnss),outStruc.plFull(idx,indsGnss),'k.','linewidth',2);
            
            ylabel(labels{idx})
            grid on
            %         ylim([0 5])
            if idx == 1
                %                 title('Position Error in ENU')
            end
        end
        legend('Error','Protection Level')
        xlabel('Time of run [min]')
    end
    
    
    % Convert position error to body frame
    if plot_pos_err_body_frame || plot_all
        figure; clf; hold on;
        errBody = nan(size(errEnu));
        for idx = 1:max(size(errBody))
            Rbodyi = euler2dcm123(truthAttEnu(idx,:));
            
            errBody(idx,:) = Rbodyi'*-errEnu(idx,:)';
            
        end
        plot(errBody)
        legend('back','right','up')
    end
    
    if plot_ground_track || plot_all
        fig2 = figure; clf; hold on;
        plot((truthEnuInterp(:,1)-truthEnuInterp(1,1))/1000,(truthEnuInterp(:,2)-truthEnuInterp(1,2))/1000,'go')
        plot((posOutEnu(:,1)-truthEnuInterp(1,1))/1000,(posOutEnu(:,2)-truthEnuInterp(1,2))/1000,'.')
        
        distanceTraveledTruth = sum(sqrt(sum(diff(truthEnuInterp).^2,2)));
        distanceTraveledEst = sum(sqrt(sum(diff(posOutEnu).^2,2)));
        
        %             disp(distanceTraveledTruth)
        %             disp(distanceTraveledEst)
    end
    
    if plot_att || plot_all
        fig3 = figure; clf; hold on;
        plot(epochsPlot,attEnu(:,1)*180/pi,'b')
        plot(epochsPlot,attEnu(:,2)*180/pi,'g')
        plot(epochsPlot,attEnu(:,3)*180/pi,'r')
        
        velMag = sqrt(sum(truthVelEnu.^2,2));
        
        truthAttEnu(repmat(velMag',1,3) < 1) = nan;
        
        plot(epochsPlot,truthAttEnu(:,1)*180/pi,'b--')
        plot(epochsPlot,truthAttEnu(:,2)*180/pi,'g--')
        plot(epochsPlot,truthAttEnu(:,3)*180/pi,'r--')
        plot(epochsPlot,velMag*10,'k','linewidth',2)
    end
    
    if plot_imu_bias || plot_all
        fig4 = figure; clf; hold on;
        plot(epochsPlot,biasSave-nanmedian(biasSave)*0)
    end
    
    if plot_vel_enu || plot_all
        fig5= figure; clf; hold on;
        plot(epochsPlot(indsGnss),errVelEnu(indsGnss,:))
        grid on;
        title('Velocity error in ENU')
    end
    
    
    
    if plot_google_maps || plot_all
        llh0 = xyz2llh(posMeas.pos)';
        llhEst = xyz2llh(pos)';
        
        fig6 = figure;; clf;
        hold on;
        plot(llh0(2,1),llh0(1,1),'go','linewidth',2)
        plot(llh0(2,end),llh0(1,end),'ro','linewidth',2)
        plot(llh0(2,:),llh0(1,:),'linewidth',2)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         plot_google_map('apiKey', 'AIzaSyAWTuImB0TjZWslTdAFiv7b5JLrPBeghFg');
        %                 plot_google_map('apiKey', 'AIzaSyAUqexHvH2M8TKgWOXwkwK4yUJHSb6wZig');
        
        %         plot_google_map;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        llhCorners = [ylim; xlim; 0 0]';
        xyzCorners = llh2xyz(llhCorners);
        
        utmZone =findUtmZone(llhCorners(1,1),llhCorners(1,2));
        
        [e1,n1,u1] = cart2utm(xyzCorners(1,1),xyzCorners(1,2),xyzCorners(1,3),utmZone);
        [e2,n2,u2] = cart2utm(xyzCorners(2,1),xyzCorners(2,2),xyzCorners(2,3),utmZone);
        
        dtick = 100;
        dticklabel = 1000;
        de0 = e2-e1;
        de = floor(de0/dtick/6)*dtick;
        dn0 = n2-n1;
        dn = floor(dn0/dtick/6)*dtick;
        
        eticks = [-3:3]*de/dticklabel;
        lonticks = [-3:3]*(de/de0)*abs(llhCorners(1,2)-llhCorners(2,2))+mean([llhCorners(1,2) llhCorners(2,2)]);
        nticks = [-3:3]*dn/dticklabel;
        latticks = [-3:3]*(dn/dn0)*abs((llhCorners(1,1)-llhCorners(2,1)))+mean([llhCorners(1,1) llhCorners(2,1)]);
        
        xticks(lonticks)
        xticklabels(num2str(eticks'))
        yticks(latticks)
        yticklabels(num2str(nticks'))
        
        xlabel('East Position [km]')
        ylabel('North Position [km]')
        
        
        drawnow;
        legend('Start','Finish')
        %         xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
    end
    
    if plot_speed || plot_all
        fig7 = figure;
        velMagKph = velMag*2.23;
        plot(epochsPlot/60,velMagKph,'k','linewidth',1);
        xlabel('Time of run [min]')
        ylabel('Vehicle speed [mph]')
    end
end

if plot_google_earth || plot_all
    
    llh0 = xyz2llh(truthXyzInterp)';
    llhEst = xyz2llh(pos)';
    
    %% build kml file
    k = kml('Pos err');
    
    f = k.createFolder('files');
    
    iconNames = repmat({''},size(llhEst,2),1);
    iconNames(1:10:size(llhEst,2)) = cellstr(num2str((1:10:size(llhEst,2))','% 4i'));
    
    alt = 0.5*ones(size(llh0(1,:)));
    
    f.plot(llh0(2,:),llh0(1,:),'altitude',alt,'linecolor','FF14F014','extrude',true);
    f.plot(llhEst(2,:),llhEst(1,:),'altitude',alt,'lineWidth',2,'linecolor','FFF02814');
    f.point(llh0(2,:),llh0(1,:),alt,'iconScale',0.5,'name','','iconColor','FF14F014');
    f.point(llhEst(2,:),llhEst(1,:),alt,'iconScale',0.5,'name',iconNames,'iconColor','FFF02814');
    
    k.run;
end


if 1
    indsGnss = find(epochSave == floor(epochSave));
end

posTrue = truthXyzInterp(indsGnss,:)';
posEst  = pos(indsGnss,:)';
epochs = epochSave(indsGnss)';
pl      = outStruc.plFull(:,indsGnss);
if isempty(outStruc.plLocFull)
    plLoc = [];
else
    plLoc   = outStruc.plLocFull(:,indsGnss);
end
errEnu  = errEnu(indsGnss,:)';





















end