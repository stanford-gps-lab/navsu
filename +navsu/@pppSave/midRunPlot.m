function midRunPlot(obj)

%% Plot covariance
figi = obj.figCov;
if ~isvalid(figi)
    obj.figResids = figure;
    figi = obj.figResids;
end
figi.Visible = 'on';
figi.Position = [figi.Position(1) figi.Position(2)+(figi.Position(4)-688) 932         688];
clf;

% List of covariance things to plot
covPlotList = {'ATTITUDE','VEL','POS','ACC_BIAS','W_BIAS','CLOCK_BIAS','TROP'};

tPlot = obj.epochSave-min(obj.epochSave);

ha = tight_subplot(3,3,0.05);
for idx = 1:length(covPlotList)
    axes(ha(idx));
    
    indsi = obj.INDS_STATE.(covPlotList{idx});
    if ~isempty(indsi)
        semilogy(tPlot,real(sqrt(obj.covSave(:,indsi))))
    end
    title(covPlotList{idx})
end


%% Plot residuals and measurement usage?

end