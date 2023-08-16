%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("complete_comparison_rounded_v2/rounded_eval_%d_1.mat",i);
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

    x_dat(i-39,1) = SNRnoised - SNRinit;
end

figure;
plot(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
plot(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
plot(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
plot(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
plot(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
plot(x_dat,snr_dat(:,6),'DisplayName','SVD (\gamma=0.0)','Color','k','LineStyle','--'); hold on
plot(x_dat,snr_dat(:,7),'DisplayName','SVD (\gamma=0.5)','Color','k','LineStyle','-'); hold on
plot(x_dat,snr_dat(:,8),'DisplayName','SVD (\gamma=1.0)','Color','k','LineStyle','-.'); hold off
legend

xlabel("Added Noise (dB)")
ylabel("SNR (dB)")

%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("complete_comparison_rounded/rounded_eval_%d_1.mat",i);
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

    x_dat(i-39,1) = SNRnoised - SNRinit;
end

figure;
plot(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
plot(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
plot(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
plot(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
plot(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
plot(x_dat,snr_dat(:,6),'DisplayName','SVD (\gamma=0.0)','Color','k','LineStyle','--'); hold on
plot(x_dat,snr_dat(:,7),'DisplayName','SVD (\gamma=0.5)','Color','k','LineStyle','-'); hold on
plot(x_dat,snr_dat(:,8),'DisplayName','SVD (\gamma=1.0)','Color','k','LineStyle','-.'); hold off
legend

xlabel("Added Noise (dB)")
ylabel("CNR (dB)")

%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("complete_comparison_rounded/rounded_eval_%d_1.mat",i);
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

    x_dat(i-39,1) = SNRnoised - SNRinit;
end

figure;
plot(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
plot(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
plot(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
plot(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
plot(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
plot(x_dat,snr_dat(:,6),'DisplayName','SVD (\gamma=0.0)','Color','k','LineStyle','--'); hold on
plot(x_dat,snr_dat(:,7),'DisplayName','SVD (\gamma=0.5)','Color','k','LineStyle','-'); hold on
plot(x_dat,snr_dat(:,8),'DisplayName','SVD (\gamma=1.0)','Color','k','LineStyle','-.'); hold off
legend

xlabel("Added Noise (dB)")
ylabel("SSIM [0,1]")

%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("complete_comparison_rounded/rounded_eval_%d_1.mat",i);
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

    x_dat(i-39,1) = SNRnoised - SNRinit;
end

figure;
plot(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
plot(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
plot(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
plot(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
plot(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
plot(x_dat,snr_dat(:,6),'DisplayName','SVD (\gamma=0.0)','Color','k','LineStyle','--'); hold on
plot(x_dat,snr_dat(:,7),'DisplayName','SVD (\gamma=0.5)','Color','k','LineStyle','-'); hold on
plot(x_dat,snr_dat(:,8),'DisplayName','SVD (\gamma=1.0)','Color','k','LineStyle','-.'); hold off
legend

xlabel("Added Noise (dB)")
ylabel("PSNR (dB)")