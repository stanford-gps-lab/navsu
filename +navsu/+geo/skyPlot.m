function skyPlot(az,el,varargin)
% Draw a skyplot given azimuth and elevation in degrees

parser = inputParser;

parser.addParameter('linestyle', '');
parser.addParameter('linecolor',[]);
parser.addParameter('label',true);
parser.addParameter('prn',[]);
parser.addParameter('const',[]);
parser.parse(varargin{:});
res = parser.Results;

linestyle = res.linestyle;
linecolor = res.linecolor;
label = res.label;
prn = res.prn;
const = res.const;

nLines = size(el,2);
for idx = 1:nLines
    elSpherical = 90*cos(el(:,idx)* pi/180);
    if isempty(linecolor)
        s = polarplot(az(:,idx)*pi/180,elSpherical*pi/180,linestyle,'linewidth',2);
    else
        s = polarplot(az(:,idx)*pi/180,elSpherical*pi/180,linestyle,'linewidth',2,'color',linecolor(idx,:));
    end
    
    % change data tips to be useful
    row = dataTipTextRow('el',el(:,idx));
    s.DataTipTemplate.DataTipRows(1) = row;
    row = dataTipTextRow('az',az(:,idx));
    s.DataTipTemplate.DataTipRows(2) = row;
        
    if ~isempty(prn) 
        row = dataTipTextRow('PRN',prn(idx).*ones(size(az(:,idx))));
        s.DataTipTemplate.DataTipRows(end+1) = row;
    end
    
    if ~isempty(const)
        row = dataTipTextRow('const',const(idx).*ones(size(az(:,idx))));
        s.DataTipTemplate.DataTipRows(end+1) = row;
    end
    
     if label && ~all(isnan(elSpherical))
        % add some text
        indi = find(~isnan(elSpherical));
        indi = indi(1);
        
        if ~isempty(prn)
            texti = num2str(prn(idx));
        else
            texti = num2str(idx);
        end
        if ~isempty(prn) && ~isempty(const)
            texti = [texti ', ' num2str(const(idx))];
        end
        
        text(az(indi,idx)*pi/180,elSpherical(indi)*pi/180,texti,'FontSize',12,...
            'FontWeight','bold')
        
     end
    hold on;
end
hAxis = gca;
hAxis.ThetaDir = 'clockwise';
hAxis.ThetaZeroLocation = 'top';

end