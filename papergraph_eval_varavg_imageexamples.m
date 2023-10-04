load("var_avg_data\compiled_result.mat");

%% plot avg 1
% cmode
im = cmodefunc(result_avg_1.SVDgamma0);
maxvalcmode = max(im,[],"all");
minvalcmode = min(im,[],"all");
%% bmode
im = bmodefunc(result_avg_1.SVDgamma0);
maxvalbmode = max(im,[],"all");
minvalbmode = min(im,[],"all");


%% plot without
titlelist = [['Before' newline 'Denoising'],"BPF (M_{raw})","SVD \gamma=0.0", "SVD \gamma=0.5","SVD \gamma=1.0","EMD-MI",...
    "Sym6-SURE","Sym4-DWT","Sym4-MODWT"];
figure("Position",[68,258,1832,643])
tl = tiledlayout(3,9,'TileSpacing','tight');

seldat = result_avg_1;
varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    line([1 100],[55 55],'LineWidth',1,'Color','g','LineStyle','--')
    line([56 76],[10 10],'LineWidth',1,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+9)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:100)*200 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+18)
taxis = (1:1024)*1e9/Fs;
dat3d = seldat.(varnames{i});

plot(taxis,dat3d(:,9,55),'DisplayName','Signal');hold on;plot(taxis,dat3d(:,37,27),'DisplayName','Noise');hold off
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([-5e-3 10e-3]);
set(gca,'FontName','Times New Roman','FontSize',12)
end


% cbh = colorbar;
% cbh.Layout.Tile = 'north';
% put legend
leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
leg.Layout.Tile = 'south';

title(tl,"Without signal averaging",'FontName','Times New Roman','FontSize',18,'FontWeight','bold')
ax = gcf;
exportgraphics(ax,'./exported_images/fig8_ex1.emf')


%% plot with 100
titlelist = [['Before' newline 'Denoising'],"BPF (M_{raw})","SVD \gamma=0.0", "SVD \gamma=0.5","SVD \gamma=1.0","EMD-MI",...
    "Sym6-SURE","Sym4-DWT","Sym4-MODWT"];
figure("Position",[68,258,1832,643])
tl = tiledlayout(3,9,'TileSpacing','tight');

seldat = result_avg_100;
varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    line([1 100],[55 55],'LineWidth',1,'Color','g','LineStyle','--')
    line([56 76],[10 10],'LineWidth',1,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+9)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:100)*200 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+18)
taxis = (1:1024)*1e9/Fs;
dat3d = seldat.(varnames{i});

plot(taxis,dat3d(:,9,55),'DisplayName','Signal');hold on;plot(taxis,dat3d(:,37,27),'DisplayName','Noise');hold off
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([-5e-3 10e-3]);
set(gca,'FontName','Times New Roman','FontSize',12)
end


% cbh = colorbar;
% cbh.Layout.Tile = 'north';
% put legend
leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
leg.Layout.Tile = 'south';

title(tl,"With 100-signal averaging",'FontName','Times New Roman','FontSize',18,'FontWeight','bold')
ax = gcf;
exportgraphics(ax,'./exported_images/fig8_ex100.emf')

%%
function p = cmodefunc(x)
    p = squeeze(max(abs(hilbert(x))));
    p = (p - min(p(:))) / range(p(:));
end

function p = bmodefunc(x)
    hb = abs(hilbert(x));
    p = squeeze(hb(:,:,55));
    p = (p - min(p(:))) / range(p(:));
end