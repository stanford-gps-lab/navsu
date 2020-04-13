function skyPlot(az,el,varargin)
% Draw a skyplot given azimuth and elevation in degrees

parser = inputParser;

parser.addParameter('linestyle', '');
parser.addParameter('linecolor',[]);
parser.addParameter('label',true);
parser.parse(varargin{:});
res = parser.Results;

linestyle = res.linestyle;
linecolor = res.linecolor;
label = res.label;

nLines = size(el,2);
for idx = 1:nLines
    elSpherical = 90*cos(el(:,idx)* pi/180);
    if isempty(linecolor)
        polarplot(az(:,idx)*pi/180,elSpherical*pi/180,linestyle,'linewidth',2)
    else
        polarplot(az(:,idx)*pi/180,elSpherical*pi/180,linestyle,'linewidth',2,'color',linecolor(idx,:))
    end
    
     if label && ~all(isnan(elSpherical))
        % add some text
        indi = find(~isnan(elSpherical));
        indi = indi(1);
        
        text(az(indi,idx)*pi/180,elSpherical(indi)*pi/180,num2str(idx))
        
     end
    hold on;
end
hAxis = gca;
hAxis.ThetaDir = 'clockwise';
hAxis.ThetaZeroLocation = 'top';

end