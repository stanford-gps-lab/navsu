function [solName,posXyz, epochs, cov, ztd, vpl,velEnu,stdEnu] = parsePosSolFile(filename)

solName = [];
posXyz  = [];
epochs  = [];
cov     = [];
ztd     = [];
vpl     = [];
velEnu  = [];
stdEnu  = [];

if ischar(filename)
    [~,~,ext] = fileparts(filename);
else
    ext = 'xyz';
end
switch ext
    case '.mat' % su ppp output
        %             solNames{idx} = 'Stanford';
        temp = load(filename);
        epochs = temp.epochsd';
        posXyz = temp.posOut;
        cov = temp.debugOutput.covPos;
        ztd  = temp.debugOutput.tropFullSave';
        posXyz(posXyz == 0) = nan;
        
        if isfield(temp.debugOutput,'resPr')
            resids.resPr = temp.debugOutput.resPr;
            resids.resPh = temp.debugOutput.resPh;
            resids.epochs = epochs;
            resids.prns  = temp.debugOutput.prns;
            resids.constInds = temp.debugOutput.constInds;
            resids.measInfoFull = temp.debugOutput.measInfoFull;
        end
        
        nBanks = size(posXyz,3);
        
        if isfield(temp.debugOutput,'vplRt')
            vpl = temp.debugOutput.vplRt;
            vpl = temp.debugOutput.pl_ss_approx1_optfa;
            
            vplRt = temp.debugOutput.vplRt;
            epochsVpls = temp.epochsd;
            vplLoc = temp.debugOutput.vplLoc;
            %                 vplLoc = temp.debugOutput.vplRt;
            pl_ss_exact_optfa   = temp.debugOutput.pl_ss_exact_optfa;
            pl_ss_approx1_optfa = temp.debugOutput.pl_ss_approx1_optfa;
            pl_ss_approx2_optfa = temp.debugOutput.pl_ss_approx2_optfa;
            pl_ss_exact         = temp.debugOutput.pl_ss_exact;
            pl_ss_approx1       = temp.debugOutput.pl_ss_approx1;
            pl_ss_approx2       = temp.debugOutput.pl_ss_approx2;
            pl_chi2_exact       = temp.debugOutput.pl_chi2_exact;
            pl_chi2_approx1     = temp.debugOutput.pl_chi2_approx1;
            pl_chi2_approx2     = temp.debugOutput.pl_chi2_approx2;
        end
        
    case '.txt'
        % novatel truth
        solName =  'Novatel Truth';
        novTruth = csvread(filename,1);
        
        epochs = navsu.time.gps2epochs(novTruth(:,1),novTruth(:,2));
        llhNov    = novTruth(:,3:5);
        posXyz = navsu.geo.llh2xyz(llhNov,1);
        velEnu = novTruth(:,10:12);
        
        cov   = zeros(3,3,length(epochs));
        cov(1,1,:) = novTruth(:,17);
        cov(2,2,:) = novTruth(:,18);
        cov(3,3,:) = novTruth(:,19);
        cov(2,1,:) = novTruth(:,20);
        cov(1,2,:) = novTruth(:,20);
        cov(3,1,:) = novTruth(:,21);
        cov(1,3,:) = novTruth(:,21);
        cov(3,2,:) = novTruth(:,22);
        cov(2,3,:) = novTruth(:,22);
        
    case '.truth'
        % Different novatel truth
        solName =  'Novatel Truth';
        novTruth = csvread(filename,1);
        
        novTruth = [2075*ones(size(novTruth,1),1) novTruth];
        
        epochs = navsu.time.gps2epochs(novTruth(:,1),novTruth(:,2));
        llhNov    = novTruth(:,3:5);
        posXyz = navsu.geo.llh2xyz(llhNov,1);
        if size(novTruth,2) > 11
            velEnu = novTruth(:,10:12);
            
            cov   = zeros(3,3,length(epochs));
            cov(1,1,:) = novTruth(:,17);
            cov(2,2,:) = novTruth(:,18);
            cov(3,3,:) = novTruth(:,19);
            cov(2,1,:) = novTruth(:,20);
            cov(1,2,:) = novTruth(:,20);
            cov(3,1,:) = novTruth(:,21);
            cov(1,3,:) = novTruth(:,21);
            cov(3,2,:) = novTruth(:,22);
            cov(2,3,:) = novTruth(:,22);
        end
        
    case '.pos2' % rtklib
        solName = 'RTK';
        
        datai = importdata(filename,' ',26);
        
        timeCells = datai.textdata(27:end,1:2);
        timeData = strcat(timeCells(:,1),repmat({' '},size(timeCells,1),1),timeCells(:,2));
        
        epochs = round(1000*datenum2epochs(datenum(timeData)))/1000;
        xyzPos = datai.data(:,1:3);
        cov = datai.data(:,6:8);
        
        % mask to only fix data
        solType = datai.data(:,4);
        
        indsKeep = find(solType == 1 | solType == 200);
        
        epochs = epochs(indsKeep);
        xyzPos =  xyzPos(indsKeep,:);
        cov = cov(indsKeep,:);
        
    case '.pos'% nrc ppp
        solName = 'NRC';
        dataNRC    = importdata(filename,' ',8);
        datestrNRC = strcat(dataNRC.textdata(9:end,5),repmat({' '},length(dataNRC.textdata(9:end,1)),1), dataNRC.textdata(9:end,6));
        
        epochs = navsu.time.datenum2epochs(datenum(datestrNRC));
        epochs = round(epochs*20)/20;
        
        % old and busted
        try
            latNRC = [cellfun(@str2num,dataNRC.textdata(9:end,21))+cellfun(@str2num,dataNRC.textdata(9:end,22))/60+cellfun(@str2num,dataNRC.textdata(9:end,23))/3600];
            lonNRC = [cellfun(@str2num,dataNRC.textdata(9:end,24))-cellfun(@str2num,dataNRC.textdata(9:end,25))/60-cellfun(@str2num,dataNRC.textdata(9:end,26))/3600];
            hNRC   = [cellfun(@str2num,dataNRC.textdata(9:end,27))];
            ztd    = cellfun(@str2num,dataNRC.textdata(9:end,15));
        catch
            % new hotness
            latNRC = [dataNRC.data(:,15)+dataNRC.data(:,16)/60+dataNRC.data(:,17)/3600];
            lonNRC = [dataNRC.data(:,18)+sign(dataNRC.data(:,18)).*dataNRC.data(:,19)/60+sign(dataNRC.data(:,18)).*dataNRC.data(:,20)/3600];
            hNRC   = [dataNRC.data(:,21)];
            ztd = dataNRC.data(:,9);
            
            stdEnu = [dataNRC.data(:,11) dataNRC.data(:,10) dataNRC.data(:,12)]/2;
        end
        llhNrc = [latNRC lonNRC hNRC];
        
        dpos = -[0.8756
            -1.4496
            0.0382];
        posXyz = navsu.geo.llh2xyz(llhNrc,1)+0*dpos';
        
        % sort everything
        datai = sortrows([epochs posXyz ztd],1);
        epochs = datai(:,1);
        posXyz = datai(:,2:4);
        ztd    = datai(:,5);
        
        %         solData = [solData; {posNrc' epochs' nan(3,3,length(epochs)) ztdNrc0' nan(size(posNrc'))}];
        
    case  '.sum' % jpl
        solNames = [solNames; 'JPL'];
        dataNRC    = importdata(filename,' ',18);
        datestrNRC = [dataNRC.textdata(19:end,2)];
        
        datestrNRC = cellfun(@(c)[c(1:4) '/' c(6:7) '/' c(9:10) ' ' c(12:end)],datestrNRC,'uni',false);
        
        epochs = datenum2epochs(datenum(datestrNRC));
        epochs = round(epochs*20)/20;
        
        posNrc = [dataNRC.data(:,1) dataNRC.data(:,3) dataNRC.data(:,5)];
        
        ztd = zeros(size(posNrc,1),1);
        
        solData = [solData; {posNrc' epochs' nan(3,3,length(epochs)) ztd' nan(size(posNrc'))}];
    case 'xyz' % IGS XYZ position
        
        if size(filename,1) > size(filename,2)
            posXyz = filename';
        else
            posXyz = filename;
        end
        
        
end

end