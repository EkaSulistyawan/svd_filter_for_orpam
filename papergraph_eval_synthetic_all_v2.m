figure('Position',[1921,123,1918,963])
Tbig = tiledlayout(2,1);

%% for flat
Tsmall = tiledlayout(Tbig,1,5);
Tsmall.Layout.Tile = 1;
%% create graph CNR

snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    %filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.snr;
    snr_dat(i-39,2) = evalPaper1EMDMI.snr;
    snr_dat(i-39,3) = evalPaper1Wavelet.snr;
    snr_dat(i-39,4) = evalPaper2DWT.snr;
    snr_dat(i-39,5) = evalPaper2MODWT.snr;
    snr_dat(i-39,6) = evalSVDgamma0.snr;
    snr_dat(i-39,7) = evalSVDgammahalf.snr;
    snr_dat(i-39,8) = evalSVDgamma1.snr;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
%legend('Location','southeast')
axis tight
ylim([0 40])

xlabel("SNR before denoising (dB)")
ylabel("SNR after denoising (dB)")

set(gca,'FontName','Times New Roman','FontSize',18)

%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.cnr;
    snr_dat(i-39,2) = evalPaper1EMDMI.cnr;
    snr_dat(i-39,3) = evalPaper1Wavelet.cnr;
    snr_dat(i-39,4) = evalPaper2DWT.cnr;
    snr_dat(i-39,5) = evalPaper2MODWT.cnr;
    snr_dat(i-39,6) = evalSVDgamma0.cnr;
    snr_dat(i-39,7) = evalSVDgammahalf.cnr;
    snr_dat(i-39,8) = evalSVDgamma1.cnr;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([-40 30])

xlabel("SNR before denoising (dB)")
ylabel("CNR after denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',18)
%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.ssim_cmode;
    snr_dat(i-39,2) = evalPaper1EMDMI.ssim_cmode;
    snr_dat(i-39,3) = evalPaper1Wavelet.ssim_cmode;
    snr_dat(i-39,4) = evalPaper2DWT.ssim_cmode;
    snr_dat(i-39,5) = evalPaper2MODWT.ssim_cmode;
    snr_dat(i-39,6) = evalSVDgamma0.ssim_cmode;
    snr_dat(i-39,7) = evalSVDgammahalf.ssim_cmode;
    snr_dat(i-39,8) = evalSVDgamma1.ssim_cmode;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([0 1])

xlabel("SNR before denoising (dB)")
ylabel("SSIM after denoising [0,1]")
set(gca,'FontName','Times New Roman','FontSize',18)
%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.psnr_cmode;
    snr_dat(i-39,2) = evalPaper1EMDMI.psnr_cmode;
    snr_dat(i-39,3) = evalPaper1Wavelet.psnr_cmode;
    snr_dat(i-39,4) = evalPaper2DWT.psnr_cmode;
    snr_dat(i-39,5) = evalPaper2MODWT.psnr_cmode;
    snr_dat(i-39,6) = evalSVDgamma0.psnr_cmode;
    snr_dat(i-39,7) = evalSVDgammahalf.psnr_cmode;
    snr_dat(i-39,8) = evalSVDgamma1.psnr_cmode;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([0 45])

xlabel("SNR before denoising (dB)")
ylabel("PSNR after denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',18)

%% create graph time
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    %filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.time;
    snr_dat(i-39,2) = evalPaper1EMDMI.time;
    snr_dat(i-39,3) = evalPaper1Wavelet.time;
    snr_dat(i-39,4) = evalPaper2DWT.time;
    snr_dat(i-39,5) = evalPaper2MODWT.time;
    snr_dat(i-39,6) = evalSVDgamma0.time;
    snr_dat(i-39,7) = evalSVDgammahalf.time;
    snr_dat(i-39,8) = evalSVDgamma1.time;

    x_dat(i-39,1) = SNRnoised;
end
ax = nexttile(Tsmall);
boxchart(snr_dat)
ylabel("Time (s)")
ylim([0 8])
set(gca,'FontName','Times New Roman','FontSize',18)
ax.XTickLabel = {'BPF','EMD-MI','Sym6-SURE','Sym4-DWT','Sym4-MODWT','SVD (\gamma=0.0)','SVD (\gamma=0.5)','SVD (\gamma=1.0)'};
ax.XAxis.FontSize = 14;

% for EMD-MI
A = mean(snr_dat);
B = std(snr_dat);
fprintf("%.2f (%.2f)\n",A(2),B(2))

title(Tsmall,"Evaluation of USAF-1951","FontName","Times New Roman","FontSize",18,"FontWeight","bold");


%% evaluation of spherical
%% for flat
Tsmall = tiledlayout(Tbig,1,5);
Tsmall.Layout.Tile = 2;

%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_rounded_woBPF/rounded_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.snr;
    snr_dat(i-39,2) = evalPaper1EMDMI.snr;
    snr_dat(i-39,3) = evalPaper1Wavelet.snr;
    snr_dat(i-39,4) = evalPaper2DWT.snr;
    snr_dat(i-39,5) = evalPaper2MODWT.snr;
    snr_dat(i-39,6) = evalSVDgamma0.snr;
    snr_dat(i-39,7) = evalSVDgammahalf.snr;
    snr_dat(i-39,8) = evalSVDgamma1.snr;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
%legend('Location','southeast')
axis tight
ylim([0 40])

xlabel("SNR before denoising (dB)")
ylabel("SNR after denoising (dB)")

set(gca,'FontName','Times New Roman','FontSize',18)

%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_rounded_woBPF/rounded_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.cnr;
    snr_dat(i-39,2) = evalPaper1EMDMI.cnr;
    snr_dat(i-39,3) = evalPaper1Wavelet.cnr;
    snr_dat(i-39,4) = evalPaper2DWT.cnr;
    snr_dat(i-39,5) = evalPaper2MODWT.cnr;
    snr_dat(i-39,6) = evalSVDgamma0.cnr;
    snr_dat(i-39,7) = evalSVDgammahalf.cnr;
    snr_dat(i-39,8) = evalSVDgamma1.cnr;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([-40 30])

xlabel("SNR before denoising (dB)")
ylabel("CNR after denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',18)
%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_rounded_woBPF/rounded_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.ssim_cmode;
    snr_dat(i-39,2) = evalPaper1EMDMI.ssim_cmode;
    snr_dat(i-39,3) = evalPaper1Wavelet.ssim_cmode;
    snr_dat(i-39,4) = evalPaper2DWT.ssim_cmode;
    snr_dat(i-39,5) = evalPaper2MODWT.ssim_cmode;
    snr_dat(i-39,6) = evalSVDgamma0.ssim_cmode;
    snr_dat(i-39,7) = evalSVDgammahalf.ssim_cmode;
    snr_dat(i-39,8) = evalSVDgamma1.ssim_cmode;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([0 1])

xlabel("SNR before denoising (dB)")
ylabel("SSIM after denoising [0,1]")
set(gca,'FontName','Times New Roman','FontSize',18)
%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_rounded_woBPF/rounded_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.psnr_cmode;
    snr_dat(i-39,2) = evalPaper1EMDMI.psnr_cmode;
    snr_dat(i-39,3) = evalPaper1Wavelet.psnr_cmode;
    snr_dat(i-39,4) = evalPaper2DWT.psnr_cmode;
    snr_dat(i-39,5) = evalPaper2MODWT.psnr_cmode;
    snr_dat(i-39,6) = evalSVDgamma0.psnr_cmode;
    snr_dat(i-39,7) = evalSVDgammahalf.psnr_cmode;
    snr_dat(i-39,8) = evalSVDgamma1.psnr_cmode;

    x_dat(i-39,1) = SNRnoised;
end

ax = nexttile(Tsmall);;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight
ylim([0 45])

xlabel("SNR before denoising (dB)")
ylabel("PSNR after denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',18)
%% legend
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'south';

set(gca,'FontName','Times New Roman','FontSize',18)
%% create graph time
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_rounded_woBPF/rounded_eval_%d_1.mat",i);
    %filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_flat/flat_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF.time;
    snr_dat(i-39,2) = evalPaper1EMDMI.time;
    snr_dat(i-39,3) = evalPaper1Wavelet.time;
    snr_dat(i-39,4) = evalPaper2DWT.time;
    snr_dat(i-39,5) = evalPaper2MODWT.time;
    snr_dat(i-39,6) = evalSVDgamma0.time;
    snr_dat(i-39,7) = evalSVDgammahalf.time;
    snr_dat(i-39,8) = evalSVDgamma1.time;

    x_dat(i-39,1) = SNRnoised;
end
ax = nexttile(Tsmall);
boxchart(snr_dat)
ylabel("Time (s)")
ylim([0 8])
set(gca,'FontName','Times New Roman','FontSize',18)
ax.XTickLabel = {'BPF','EMD-MI','Sym6-SURE','Sym4-DWT','Sym4-MODWT','SVD (\gamma=0.0)','SVD (\gamma=0.5)','SVD (\gamma=1.0)'};
ax.XAxis.FontSize = 14;

% for EMD-MI
A = mean(snr_dat);
B = std(snr_dat);
fprintf("%.2f (%.2f)\n",A(2),B(2))
title(Tsmall,"Evaluation of Red-stained Spherical Phantom","FontName","Times New Roman","FontSize",18,"FontWeight","bold");


exportgraphics(gcf,'./exported_images/fig7_eval_synthetic.emf')