function obsOut = stripMeas(index,obsFull,obsMap)

obsMapi = obsMap(index,:);
obsOut = [];

for idx = 1:length(obsMapi)
    if obsMapi(idx) ~= 0
        obsFulli = obsFull{idx};
        
        switch obsFulli.type
            case navsu.internal.MeasEnum.GNSS
                obsi = obsFulli;
                obsi.epochs = obsFulli.epochs(obsMapi(idx));
                obsi.range.obs   = obsFulli.range.obs(:,:,obsMapi(idx));
                obsi.doppler.obs = obsFulli.doppler.obs(:,:,obsMapi(idx));
                obsi.snr.obs     = obsFulli.snr.obs(:,:,obsMapi(idx));                
                
                % Receiver output time of lock
                if ~isempty(obsi.range.lockTime)
                    obsi.range.lockTime = obsi.range.lockTime(:,:,obsMapi(idx));
                end
                
                
                
            case navsu.internal.MeasEnum.Position
                obsi = [];
                obsi.epochs = obsFulli.epochs(obsMapi(idx));
                obsi.obs    = obsFulli.obs(:,obsMapi(idx));
                obsi.cov    = obsFulli.cov(:,:,obsMapi(idx));
                obsi.type   = obsFulli.type;
                obsi.ID     = obsFulli.ID;
                
            case navsu.internal.MeasEnum.Velocity
                obsi = [];
                obsi.epochs = obsFulli.epochs(obsMapi(idx));
                obsi.obs    = obsFulli.obs(:,obsMapi(idx));
                obsi.cov    = obsFulli.cov(:,:,obsMapi(idx));
                obsi.type   = obsFulli.type;
                obsi.ID     = obsFulli.ID;
                
            case navsu.internal.MeasEnum.IMU
                imuFields = fields(obsFulli);
                for fdx = 1:length(imuFields)
                    if size(obsFulli.(imuFields{fdx}),1) > 1
                        obsi.(imuFields{fdx}) = obsFulli.(imuFields{fdx})(obsMapi(idx),:);
                    else
                        obsi.(imuFields{fdx}) = obsFulli.(imuFields{fdx});
                    end
                end
                
            case navsu.internal.MeasEnum.Wheels
                imuFields = fields(obsFulli);
                for fdx = 1:length(imuFields)
                    if size(obsFulli.(imuFields{fdx}),1) > 1
                        obsi.(imuFields{fdx}) = obsFulli.(imuFields{fdx})(obsMapi(idx),:);
                    else
                        obsi.(imuFields{fdx}) = obsFulli.(imuFields{fdx});
                    end
                end
                
        end
        obsOut = [obsOut; {obsi}];
    end
end


end