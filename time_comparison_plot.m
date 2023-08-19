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
boxchart(snr_dat)
set(gca,'XTickLabel',{'BPF','EMD-MI','Sym6-SURE','Sym4-DWT','Sym4-MODWT','SVD (\gamma=0.0)','SVD (\gamma=0.5)','SVD (\gamma=1.0)'});
ylabel("Time (s)")
ylim([0 8])
set(gca,'FontName','Times New Roman','FontSize',14)
% scatter(x_dat,snr_dat(:,1),'DisplayName','BPF'); hold on
% scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI'); hold on
% scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE'); hold on
% scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT'); hold on
% scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT'); hold on
% scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
% scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
% scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
% %legend('Location','southeast')