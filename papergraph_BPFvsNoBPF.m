figure("Position",[552,458,1244,420])
ax = tiledlayout(1,2);
%% create graph CNR
snr_dat = zeros(40,8);
x_dat = zeros(40,1);

for i=40:80
    filename = sprintf("From DL PC/svd_filter_for_orpam/complete_comparison_rounded_woBPF/rounded_eval_%d_1.mat",i);
    load(filename)
    % load metrics
    snr_dat(i-39,1) = evalBPF_wobpf.snr;
    snr_dat(i-39,2) = evalPaper1EMDMI_wobpf.snr;
    snr_dat(i-39,3) = evalPaper1Wavelet_wobpf.snr;
    snr_dat(i-39,4) = evalPaper2DWT_wobpf.snr;
    snr_dat(i-39,5) = evalPaper2MODWT_wobpf.snr;
    snr_dat(i-39,6) = evalSVDgamma0_wobpf.snr;
    snr_dat(i-39,7) = evalSVDgammahalf_wobpf.snr;
    snr_dat(i-39,8) = evalSVDgamma1_wobpf.snr;

    x_dat(i-39,1) = SNRnoised;
end

nexttile();
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI','MarkerEdgeColor',[0.8500 0.3250 0.0980]); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE','MarkerEdgeColor',[0.9290 0.6940 0.1250]); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT','MarkerEdgeColor',[0.4940 0.1840 0.5560]); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT','MarkerEdgeColor',[0.4660 0.6740 0.1880]); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
set(gca,"FontName","Times New Roman","FontSize",12)
title("Performance without prior BPF",'FontName','Times New Roman','FontSize',14)
xlabel("SNR prior denoising (dB)",'FontName','Times New Roman','FontSize',14)
ylabel("SNR after denoising (dB)",'FontName','Times New Roman','FontSize',14)
axis tight
ylim([-5 45])
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

nexttile();
scatter(x_dat,snr_dat(:,2),'DisplayName','EMD-MI','MarkerEdgeColor',[0.8500 0.3250 0.0980]); hold on
scatter(x_dat,snr_dat(:,3),'DisplayName','Sym6-SURE','MarkerEdgeColor',[0.9290 0.6940 0.1250]); hold on
scatter(x_dat,snr_dat(:,4),'DisplayName','Sym4-DWT','MarkerEdgeColor',[0.4940 0.1840 0.5560]); hold on
scatter(x_dat,snr_dat(:,5),'DisplayName','Sym4-MODWT','MarkerEdgeColor',[0.4660 0.6740 0.1880]); hold on
scatter(x_dat,snr_dat(:,6),[],'k','+','DisplayName','SVD (\gamma=0.0)'); hold on
scatter(x_dat,snr_dat(:,7),[],'k','o','DisplayName','SVD (\gamma=0.5)'); hold on
scatter(x_dat,snr_dat(:,8),[],'k','x','DisplayName','SVD (\gamma=1.0)'); hold off
set(gca,"FontName","Times New Roman","FontSize",12)
title("Performance with prior BPF",'FontName','Times New Roman','FontSize',14)
xlabel("SNR prior denoising (dB)",'FontName','Times New Roman','FontSize',14)
ylabel("SNR after denoising (dB)",'FontName','Times New Roman','FontSize',14)
axis tight
ylim([-5 45])
%%
leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
leg.Layout.Tile = 'south';


ax = gcf;
exportgraphics(ax,'./exported_images/fig12_compare.emf')