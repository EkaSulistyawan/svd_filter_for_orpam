load("saved_rbc\rbc2.mat");

%%
figure
imagesc(cmodefunc(space3Dbpf)')
colormap hot
axis image off tight
ax = gcf;
exportgraphics(ax,'./exported_images/rbc_bpf.emf')

figure
imagesc(cmodefunc(denoisedSVDgammahalf)')
colormap hot
axis image off tight
ax = gcf;
exportgraphics(ax,'./exported_images/rbc_gammahalf.emf')

%%
function p = cmodefunc(x)
    p = squeeze(max(abs(hilbert(x))));
    p = (p - min(p(:))) / range(p(:));
    % p(p > 0.6) = 0.6;
    % p = (p - min(p(:))) / range(p(:));
end

function p = bmodefunc(x)
    hb = abs(hilbert(x));
    p = squeeze(hb(:,:,71));
    p = (p - min(p(:))) / range(p(:));
end