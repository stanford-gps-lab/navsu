function [trop0,m,tropDataSave] = tropDelay(el,az, h, lat, lon, doy, PARAMS, tropData,INDS_STATIONS,epochi,preloadData)

% Elevation and azimuth are in degrees

switch PARAMS.tropModel
    case 'UNB3'
        [trop0,m,tropDataSave] = navsu.ppp.models.tropoErrCorrUnb3(el,h,lat,doy);
        
    case 'IGS'
        [trop0,m,tropDataSave] = tropo_error_correction_IGS(el,az,h,lat,lon,doy,PARAMS,tropData,...
            INDS_STATIONS,epochi);
        
    case 'SAAS'
        % saastomoinen
        [trop0] = saastamoinen_model_SU(lat, lon, h, el);
        m = zeros(size(trop0));
        tropDataSave.trototSave = zeros(size(trop0));
        tropDataSave.gmfwSave  = zeros(size(trop0));
        
    case 'GO'
        [trop0,m,tropDataSave] = tropo_error_correction_go(el,h,lat,doy);
        
        
    case 'GMF' 
        mjd = jd2mjd(epochs2jd(epochi));
        
        [trop0,m,tropDataSave] = tropo_error_correction_gmf(el,h,lat,lon,mjd,epochi);
end

end