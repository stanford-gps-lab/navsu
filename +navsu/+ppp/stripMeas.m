function obsi = stripMeas(obsFull,obsType,index)

obsList = obsFull{obsType};

switch obsType
    case 1
        % GNSS
        obsi = obsList;
        obsi.epochs = obsList.epochs(index);
        obsi.range.obs   = obsList.range.obs(:,:,index);
        obsi.doppler.obs = obsList.doppler.obs(:,:,index);
        obsi.snr.obs     = obsList.snr.obs(:,:,index);
        
        % Receiver output time of lock
        if ~isempty(obsi.range.lockTime)
            obsi.range.lockTime = obsi.range.lockTime(:,:,index);
        end
        
    case 2
        % IMU
        imuFields = fields(obsList);
        for idx = 1:length(imuFields)
           obsi.(imuFields{idx}) = obsList.(imuFields{idx})(index,:); 
        end
end


end