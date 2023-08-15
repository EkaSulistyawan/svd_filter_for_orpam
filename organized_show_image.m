function [outp] = organized_show_image(dat)
% for tiledLayout resources: https://www.mathworks.com/matlabcentral/answers/1444554-how-to-set-common-colorbar-for-multiplots
% before denoising
hb = abs(hilbert(dat));
maxHb = max(hb,[],"all");
minHb = min(hb,[],"all");

figure;
tcl = tiledlayout(2,1);
cmode = squeeze(max(hb));
nexttile();imagesc(cmode');axis off;colormap hot;clim([minHb maxHb])

where_slice = 14; % on USAF
where_slice = 163; % on 09102022_sphere
bmode = squeeze(hb(:,:, where_slice));
nexttile();imagesc(bmode);axis off;colormap hot;clim([minHb maxHb])

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Intensity (a.u.)';

set(gca,'FontName','Times New Roman','FontSize',15,'FontWeight','bold')
objects = get(gcf, 'Children');
currentpos = get(objects(1),'Position');
currentpos(3) = currentpos(3) - 0.1;
set(objects(1), 'Position', currentpos);
set(gcf,'Position',[400 100 403 800])
outp =1;
end

