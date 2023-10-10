%% USAF-1951 & Red-stained spherical phantom

figure('Position',[680,458,1049,420])
global_param
ax = tiledlayout(1,2);

nexttile()
load("../../SharingPoint/data/cellular/compiled/04092023_USAF.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:100,1:100);
Fs = 5e9;
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);
cmode = squeeze(max(abs(hilbert(space3D))));
imagesc(cmode');axis image off;colormap hot
%line([156 176],[188 188],'LineWidth',2,'Color','white')

global_param;
sgx = sgx_flat2; 
sgy = sgy_flat2;
bgx = bgx_flat2; 
bgy = bgy_flat2;
rectangle('Position',[sgx(1) sgy(1) range(sgx) range(sgy)],'LineWidth',2,'EdgeColor','b')
rectangle('Position',[bgx(1) bgy(1) range(bgx) range(bgy)],'LineWidth',2,'EdgeColor','g')
title("ROI for USAF-1951",'FontName','Times New Roman','FontSize',14)

nexttile()
load("../../SharingPoint/data/cellular/compiled/09102022_sphere.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:200,1:200);
Fs = 5e9;
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);
cmode = squeeze(max(abs(hilbert(space3D))));
imagesc(cmode');axis image off;colormap hot
%line([156 176],[100 100],'LineWidth',2,'Color','white')
global_param;
sgx = sgx_round; 
sgy = sgy_round;
bgx = bgx_round; 
bgy = bgy_round;
rectangle('Position',[sgx(1) sgy(1) range(sgx) range(sgy)],'LineWidth',2,'EdgeColor','b')
rectangle('Position',[bgx(1) bgy(1) range(bgx) range(bgy)],'LineWidth',2,'EdgeColor','g')
title("ROI for Red-stained Spherical Phantom",'FontName','Times New Roman','FontSize',14)

ax = gcf;
exportgraphics(ax,'./exported_images/fig7_ROI_flat2.emf')

%% fig 9 roi
load("var_avg_data\compiled_result.mat");
cmode = cmodefunc(result_avg_100.raw);
figure;
imagesc(cmode');axis image off;colormap hot
%line([156 176],[100 100],'LineWidth',2,'Color','white')
global_param;
sgx = sgx_varavg; 
sgy = sgy_varavg;
bgx = bgx_varavg; 
bgy = bgy_varavg;
rectangle('Position',[sgx(1) sgy(1) range(sgx) range(sgy)],'LineWidth',2,'EdgeColor','b')
rectangle('Position',[bgx(1) bgy(1) range(bgx) range(bgy)],'LineWidth',2,'EdgeColor','g')

ax = gcf;
exportgraphics(ax,'./exported_images/fig9_ROI.emf')

%% fig 11
load("saved_rbc\rbc.mat");

%%
cmode = cmodefunc(space3D);
figure;
imagesc(cmode');axis image off;colormap hot
%line([156 176],[100 100],'LineWidth',2,'Color','white')
%global_param;
sgx = sgx_rbc; 
sgy = sgy_rbc;
bgx = bgx_rbc; 
bgy = bgy_rbc;
rectangle('Position',[sgx(1) sgy(1) range(sgx) range(sgy)],'LineWidth',2,'EdgeColor','b')
rectangle('Position',[bgx(1) bgy(1) range(bgx) range(bgy)],'LineWidth',2,'EdgeColor','g')

ax = gcf;
exportgraphics(ax,'./exported_images/fig11_ROI.emf')

%%
cmode = cmodefunc(space3D);
figure;
imagesc(cmode');axis image off;colormap hot
%line([156 176],[100 100],'LineWidth',2,'Color','white')
global_param;
sgx = sgx_mel; 
sgy = sgy_mel;
bgx = bgx_mel; 
bgy = bgy_mel;
rectangle('Position',[sgx(1) sgy(1) range(sgx) range(sgy)],'LineWidth',2,'EdgeColor','b')
rectangle('Position',[bgx(1) bgy(1) range(bgx) range(bgy)],'LineWidth',2,'EdgeColor','g')

ax = gcf;
exportgraphics(ax,'./exported_images/melanoma_roi.emf')

%%
cmode = cmodefunc(space3D);
figure;
imagesc(cmode');axis image off;colormap hot
%line([156 176],[100 100],'LineWidth',2,'Color','white')
global_param;
sgx = sgx_rbcArt; 
sgy = sgy_rbcArt;
bgx = bgx_rbcArt; 
bgy = bgy_rbcArt;
rectangle('Position',[sgx(1) sgy(1) range(sgx) range(sgy)],'LineWidth',2,'EdgeColor','b')
rectangle('Position',[bgx(1) bgy(1) range(bgx) range(bgy)],'LineWidth',2,'EdgeColor','g')

ax = gcf;
exportgraphics(ax,'./exported_images/rbcart_roi.emf')
%% functions
function p = cmodefunc(x)
    p = squeeze(max(abs(hilbert(x))));
    p = (p - min(p(:))) / range(p(:));
    p(p > 0.45) = 0.45;
    p = (p - min(p(:))) / range(p(:));
end