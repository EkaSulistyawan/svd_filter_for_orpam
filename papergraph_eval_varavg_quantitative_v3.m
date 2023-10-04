% figure('Position',[1827,103,1640,963]); % v1
%figure('Position',[2009,143,1640,914]); % v2
figure('Position',[1921,123,1920,963]); % v3
% 
tiledlayout(2,4,'TileSpacing','loose')
%% show roi
ax = nexttile(1,[2 1])
load("var_avg_data\compiled_result.mat");
cmode = cmodefunc(result_avg_100.raw);
imagesc(cmode');axis image off;
colormap hot
cbh = colorbar;
cbh.Label.String="Normalized Intensity [0, 1]";
cbh.Location = 'northoutside';
set(gca,"Visible",false,"FontName","Times New Roman","FontSize",18)
%line([156 176],[100 100],'LineWidth',2,'Color','white')
global_param;
sgx = sgx_varavg; 
sgy = sgy_varavg;
bgx = bgx_varavg; 
bgy = bgy_varavg;
rectangle('Position',[sgx(1) sgy(1) range(sgx) range(sgy)],'LineWidth',2,'EdgeColor','b')
rectangle('Position',[bgx(1) bgy(1) range(bgx) range(bgy)],'LineWidth',2,'EdgeColor','g')

%% create evolution
x_dat = zeros(100,1);

for i=1:100
    filename = sprintf("From DL PC/svd_filter_for_orpam/comparison_varying_avg/rounded_eval_%d.mat",i);
    load(filename)

    x_dat(i,1) = SNRinit;
end

ax = nexttile;
plot([1:100],x_dat,'Color','r','LineWidth',2);
%legend('Location','southeast')
axis tight
ylim([0 inf])
title("(A)")
xlabel("Number of averaged signal")
ylabel("SNR before denoising (dB)")

set(gca,'FontName','Times New Roman','FontSize',18)

%% create graph CNR
snr_dat = zeros(100,8);
x_dat = zeros(100,1);

for i=1:100
    filename = sprintf("From DL PC/svd_filter_for_orpam/comparison_varying_avg/rounded_eval_%d.mat",i);
    load(filename)
    % load metrics
    snr_dat(i,1) = evalBPF.snr;
    snr_dat(i,2) = evalPaper1EMDMI.snr;
    snr_dat(i,3) = evalPaper1Wavelet.snr;
    snr_dat(i,4) = evalPaper2DWT.snr;
    snr_dat(i,5) = evalPaper2MODWT.snr;
    snr_dat(i,6) = evalSVDgamma0.snr;
    snr_dat(i,7) = evalSVDgammahalf.snr;
    snr_dat(i,8) = evalSVDgamma1.snr;

    x_dat(i,1) = SNRinit;
end

ax = nexttile;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF (M_{raw})'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
%legend('Location','southeast')
axis tight
ylim([0 inf])

title("(B)")
xlabel("SNR before denoising (dB)")
ylabel("SNR after denoising (dB)")

set(gca,'FontName','Times New Roman','FontSize',18)



%% create graph CNR
snr_dat = zeros(100,8);
x_dat = zeros(100,1);

for i=1:100
    filename = sprintf("From DL PC/svd_filter_for_orpam/comparison_varying_avg/rounded_eval_%d.mat",i);
    load(filename)
    % load metrics
    snr_dat(i,1) = evalBPF.cnr;
    snr_dat(i,2) = evalPaper1EMDMI.cnr;
    snr_dat(i,3) = evalPaper1Wavelet.cnr;
    snr_dat(i,4) = evalPaper2DWT.cnr;
    snr_dat(i,5) = evalPaper2MODWT.cnr;
    snr_dat(i,6) = evalSVDgamma0.cnr;
    snr_dat(i,7) = evalSVDgammahalf.cnr;
    snr_dat(i,8) = evalSVDgamma1.cnr;

    x_dat(i,1) = SNRinit;
end

ax = nexttile;;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF (M_{raw})'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([0 inf])

title("(C)")
xlabel("SNR before denoising (dB)")
ylabel("CNR after denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',18)
%% create graph CNR
snr_dat = zeros(100,8);
x_dat = zeros(100,1);

for i=1:100
    filename = sprintf("From DL PC/svd_filter_for_orpam/comparison_varying_avg/rounded_eval_%d.mat",i);
    load(filename)
    % load metrics
    snr_dat(i,1) = evalBPF.ssim_cmode;
    snr_dat(i,2) = evalPaper1EMDMI.ssim_cmode;
    snr_dat(i,3) = evalPaper1Wavelet.ssim_cmode;
    snr_dat(i,4) = evalPaper2DWT.ssim_cmode;
    snr_dat(i,5) = evalPaper2MODWT.ssim_cmode;
    snr_dat(i,6) = evalSVDgamma0.ssim_cmode;
    snr_dat(i,7) = evalSVDgammahalf.ssim_cmode;
    snr_dat(i,8) = evalSVDgamma1.ssim_cmode;

    x_dat(i,1) = SNRinit;
end

ax = nexttile;;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF (M_{raw})'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([0 inf])

title("(D)")
xlabel("SNR before denoising (dB)")
ylabel("SSIM after denoising [0,1]")
set(gca,'FontName','Times New Roman','FontSize',18)
%% create graph CNR
snr_dat = zeros(100,8);
x_dat = zeros(100,1);

for i=1:100
    filename = sprintf("From DL PC/svd_filter_for_orpam/comparison_varying_avg/rounded_eval_%d.mat",i);
    load(filename)
    % load metrics
    snr_dat(i,1) = evalBPF.psnr_cmode;
    snr_dat(i,2) = evalPaper1EMDMI.psnr_cmode;
    snr_dat(i,3) = evalPaper1Wavelet.psnr_cmode;
    snr_dat(i,4) = evalPaper2DWT.psnr_cmode;
    snr_dat(i,5) = evalPaper2MODWT.psnr_cmode;
    snr_dat(i,6) = evalSVDgamma0.psnr_cmode;
    snr_dat(i,7) = evalSVDgammahalf.psnr_cmode;
    snr_dat(i,8) = evalSVDgamma1.psnr_cmode;

    x_dat(i,1) = SNRinit;
end

ax = nexttile;;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF (M_{raw})'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([0 inf])

title("(E)")
xlabel("SNR before denoising (dB)")
ylabel("PSNR after denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',18)

%% legend
leg = legend('Orientation', 'Vertical');
leg.Layout.Tile = 'east';

set(gca,'FontName','Times New Roman','FontSize',18)

%% create graph time
snr_dat = zeros(100,8);
x_dat = zeros(100,1);

for i=1:100
    filename = sprintf("From DL PC/svd_filter_for_orpam/comparison_varying_avg/rounded_eval_%d.mat",i);
    load(filename)
    % load metrics
    snr_dat(i,1) = evalBPF.time;
    snr_dat(i,2) = evalPaper1EMDMI.time;
    snr_dat(i,3) = evalPaper1Wavelet.time;
    snr_dat(i,4) = evalPaper2DWT.time;
    snr_dat(i,5) = evalPaper2MODWT.time;
    snr_dat(i,6) = evalSVDgamma0.time;
    snr_dat(i,7) = evalSVDgammahalf.time;
    snr_dat(i,8) = evalSVDgamma1.time;

    x_dat(i,1) = SNRinit;
end

ax = nexttile;
boxchart(snr_dat)
title("(F)")
ylabel("Time (s)")
ylim([0 8])
set(gca,'FontName','Times New Roman','FontSize',18)
ax.XTickLabel = {'BPF (Mraw)','EMD-MI','Sym6-SURE','Sym4-DWT','Sym4-MODWT','SVD (\gamma=0.0)','SVD (\gamma=0.5)','SVD (\gamma=1.0)'};
ax.XAxis.FontSize = 14;

% for EMD-MI
A = mean(snr_dat);
B = std(snr_dat);
fprintf("%.2f (%.2f)\n",A(2),B(2))

exportgraphics(gcf,'./exported_images/eval_varavg_v3.emf')

function p = cmodefunc(x)
    p = squeeze(max(abs(hilbert(x))));
    p = (p - min(p(:))) / range(p(:));
    %p(p > 0.45) = 0.45;
    %p = (p - min(p(:))) / range(p(:));
end