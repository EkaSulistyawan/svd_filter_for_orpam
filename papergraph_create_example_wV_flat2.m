%% load data 
close all
clear all
load("../../SharingPoint/data/cellular/compiled/04092023_USAF.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:100,1:100);
Fs = 5e9;
%% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);


%% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,100*100);
[u,s,v] = svd(dat2d,"econ");

v3d = reshape(v,100,100,1024);
objectness = get_dctweight(v3d);

%% tile layout 1
figure("Position",[2,42,364,954]);
tiledlayout(5,2,"Padding","compact")
nexttile(1,[2 2])
imagesc(v);
set(gca,"FontName","Times New Roman","FontSize",14)

nexttile(5)
plot(v(:,1))
set(gca,"FontName","Times New Roman","FontSize",14)

nexttile(6)
plot(v(:,9))
set(gca,"FontName","Times New Roman","FontSize",14)

dearrange = reshape(v(:,1),100,100);
dctdearrange = dct2(dearrange);
nexttile(7)
imagesc(dearrange');axis off
set(gca,"FontName","Times New Roman","FontSize",14)
nexttile(9)
imagesc(dctdearrange)
set(gca,"FontName","Times New Roman","FontSize",14)

dearrange = reshape(v(:,9),100,100);
dctdearrange = dct2(dearrange);
nexttile(8)
imagesc(dearrange');axis off
set(gca,"FontName","Times New Roman","FontSize",14)
nexttile(10)
imagesc(dctdearrange)
set(gca,"FontName","Times New Roman","FontSize",14)

%% layout 2
figure("Position",[13,42,842,954]);
tiledlayout(3,2,"Padding","compact")

nexttile
maxval = max(v(:,1),[],"all");
minval = min(v(:,9),[],"all");
plot(v(:,1))
ylim([minval maxval])
ylabel("Magnitude","FontWeight","bold");
title("V Column 1")
set(gca,"FontName","Times New Roman","FontSize",14)

nexttile
plot(v(:,9))
ylim([minval maxval])
title("V Column 2")
set(gca,"FontName","Times New Roman","FontSize",14)


nexttile(3)
dearrange = reshape(v(:,1),100,100);
maxval = max(dearrange,[],"all");
minval = min(dearrange,[],"all");
imagesc(dearrange');
clim([minval maxval])
ylabel("De-arranged to 2D","FontWeight","bold");
set(gca,"FontName","Times New Roman","FontSize",14)

jet1 = nexttile(5);
dctdearrange = dct2(dearrange);
maxvaldct = max(dctdearrange,[],"all");
minvaldct = min(dctdearrange,[],"all");
imagesc(dctdearrange)
annotation("textarrow",[0.2 0.13],[0.2 0.29],'String',...
    sprintf("Low, non-flat \nDCT coefficient\n for capturing \nobject features"),...
    'FontName','Times New Roman','FontSize',13,'Color',[0 1 0])
clim([minvaldct maxvaldct])
ylabel("DCT Coefficient","FontWeight","bold");
set(gca,"FontName","Times New Roman","FontSize",14)
colormap(jet1,gray)

nexttile(4)
dearrange = reshape(v(:,9),100,100);
imagesc(dearrange');
clim([minval maxval])
cbh = colorbar;
cbh.Label.String = "Intensity (a.u.)";
set(gca,"FontName","Times New Roman","FontSize",14)

jet2 = nexttile(6);
dctdearrange = dct2(dearrange);
imagesc(mag2db(dctdearrange.^2)/2)
clim([minvaldct maxvaldct])
colormap(jet2,gray)
cbh = colorbar;
cbh.Label.String = "Intensity (a.u.)";
set(gca,"FontName","Times New Roman","FontSize",14)
exportgraphics(gcf,'./exported_images/fig6_wV_flat2.emf')