hf = figure('Units','normalized'); 
colormap hot
cbh = colorbar;
cbh.Label.String="Normalized Intensity [0, 1]";
cbh.Location = 'south';
set(gca,"Visible",false,"FontName","Times New Roman","FontSize",30)
ax = gcf;
exportgraphics(ax,'./exported_images/colorbar_size_big.emf')