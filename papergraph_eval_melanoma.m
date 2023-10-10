load("saved_others\melanoma_test.mat");


%% plot
titlelist = [['Before' newline 'Denoising'],"BPF (M_{raw})","SVD \gamma=0.0", "SVD \gamma=0.5","SVD \gamma=1.0","EMD-MI",...
    "Sym6-SURE"];
%figure("Position",[2012,113,1623,916])
figure("Position",[1977,310,1729,766])
tl = tiledlayout(3,7,'TileSpacing','tight');

seldat.raw = space3D;
seldat.bpf = space3Dbpf;
seldat.denoisedsvdgamma0 = denoisedSVDgamma0;
seldat.denoisedsvdgammahalf = denoisedSVDgammahalf;
seldat.denoisedsvdgamma1 = denoisedSVDgamma1;
seldat.denoisedpaper1emdmi = denoisedPaper1EMDMI;
seldat.denoisedpaper1wavelet = denoisedPaper1Wavelet;

varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    line([1 125],[80 80],'LineWidth',2,'Color','g','LineStyle','--')
    % line([1 33],[73 73],'LineWidth',2,'Color','c','LineStyle','--')
    line([100 120],[10 10],'LineWidth',2,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)

% nexttile(i+14)
% xaxis = (1:33)*250 / 1000;
% plot(xaxis,im(1:33,71))
% axis tight
% ylim([0.05 0.8])
% xlabel(['Length (' char(181) 'm)'])
% ylabel(' Norm. Intensity (a.u.)')
% set(gca,'FontName','Times New Roman','FontSize',12)


nexttile(i+7)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:120)*250 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+14)
taxis = (1:1024)*1e9/Fs;
dat3d = seldat.(varnames{i});
plot(taxis,dat3d(:,9,55),'DisplayName','Signal');hold on;plot(taxis,dat3d(:,37,27),'DisplayName','Noise');hold off
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([-3e-4 5e-4]);
set(gca,'FontName','Times New Roman','FontSize',12)


end


% cbh = colorbar;
% cbh.Layout.Tile = 'north';
% put legend
leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
leg.Layout.Tile = 'south';
ax = gcf;
exportgraphics(ax,'./exported_images/supp_melanoma.emf')
%%
function p = cmodefunc(x)
    p = squeeze(max(abs(hilbert(x))));
    % p = (p - min(p(:))) / range(p(:));
    % p(p > 0.45) = 0.45;
    % p = (p - min(p(:))) / range(p(:));
end

function p = bmodefunc(x)
    hb = abs(hilbert(x));
    p = squeeze(hb(:,:,80));
    p = (p - min(p(:))) / range(p(:));
end