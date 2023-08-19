figure('Position',[0,524,1916,471])
tiledlayout(1,4)
sgtitle("Evaluation on round-surfaced object",'FontName','Times New Roman','FontSize',16,'FontWeight','bold')

%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_rounded_woBPF/rounded_eval_%d_1.mat",i)
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

ax = nexttile;
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

xlabel("SNR before Denoising (dB)")
ylabel("SNR after Denoising (dB)")

set(gca,'FontName','Times New Roman','FontSize',14)

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

ax = nexttile;;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight

xlabel("SNR before Denoising (dB)")
ylabel("CNR after Denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',14)
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

ax = nexttile;;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight

xlabel("SNR before Denoising (dB)")
ylabel("SSIM after Denoising [0,1]")
set(gca,'FontName','Times New Roman','FontSize',14)
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

ax = nexttile;;
scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
axis tight

xlabel("SNR before Denoising (dB)")
ylabel("PSNR after Denoising (dB)")
set(gca,'FontName','Times New Roman','FontSize',14)

%% legend
leg = legend('Orientation', 'Horizontal');
leg.Layout.Tile = 'south';

set(gca,'FontName','Times New Roman','FontSize',14)