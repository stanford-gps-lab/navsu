sp3_fn = 'igs19362.sp3c';
if ~isfile(sp3_fn)
  url = 'https://raw.githubusercontent.com/geospace-code/georinex/refs/heads/main/src/georinex/tests/data/igs19362.sp3c';
  sp3_fn = websave(sp3_fn, url);
end

dat = navsu.readfiles.readSp3(sp3_fn);

% 3D ECEF positions
fig1 = figure(Name="SP3 ECEF Positions: " + sp3_fn);
ax1 = axes(fig1, NextPlot='add');
scatter3(ax1, dat.position(:, 1), dat.position(:, 2), dat.position(:, 3), '.')
axis(ax1, 'equal')
grid(ax1, 'on')
xlabel(ax1, 'X ECEF (m)')
ylabel(ax1, 'Y ECEF (m)')
zlabel(ax1, 'Z ECEF (m)')
title(ax1, "SP3 => ECEF Positions: " + string(dat.DateTime))

%% add translucent Earth globe for context

% create spherical Earth surface (WGS84 mean radius ~6371 km)
R = 6371000; % meters

[xs, ys, zs] = sphere(60);
xs = xs * R; ys = ys * R; zs = zs * R;

hSurf = surf(ax1, xs, ys, zs, FaceColor=[0.6 0.8 1], FaceAlpha=0.3, ...
    EdgeColor='none', DisplayName='Earth');

% improve appearance and lighting
colormap(ax1, summer)
light(ax1, 'Position',[1 0 0],'Style','infinite')
material(ax1, 'dull')
view(ax1, 3)
hold(ax1, 'off')
