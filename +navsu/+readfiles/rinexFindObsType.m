function [obsColumns,nObsTypes,obs_types_u ] = rinexFindObsType(obs_type,sysId)
% rinexFindObsTypes
% DESCRIPTION:
%   Sorry this is a bad description. This is just an internal function to
%   loadRinexObs to collect all observation types. 
% INPUT:
%   
% OUTPUT:
%   
%
% See also: navsu.readfiles.loadRinexObs

if isempty(sysId) % RINEX v2.xx
	stri = obs_type{1};

    obs_type_i = cellstr(reshape(stri,2,size(stri,2)/2)');

    obsColumns = obs_type_i;
    nObsTypes = length(obs_type_i);
    
    obs_types_u = obs_type_i;
else
    obs_types_u = {};
    for sdx = 1:length(sysId)
        stri = obs_type.(sysId{sdx});
        obs_type_i = cellstr(reshape(stri,3,size(stri,2)/3)');
        
        obsColumns.(sysId{sdx}) = obs_type_i;
        
        nObsTypes.(sysId{sdx}) = length(obs_type_i);
        
        obs_types_u = [obs_types_u; obs_type_i];
        
    end
    obs_types_u = unique(obs_types_u);
    
end
end
