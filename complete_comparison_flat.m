%% Comparison only on one signal
clearvars -except dataSNR idx ievalstart jevalstart

% USAF1
% load("../../SharingPoint/data/cellular/compiled/17122022_USAF.mat")

% USAF2
load("../../SharingPoint/data/cellular/compiled/04092023_USAF.mat")

sz = size(space3D,2);
space3D = space3D(:,1:sz-1,1:sz-1);
Fs = 5e9;

% for nonflat
global_param;
sgx = sgx_flat2; 
sgy = sgy_flat2;
bgx = bgx_flat2; 
bgy = bgy_flat2;

sg = reshape(space3D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRinit = snr2d(sg,bg);

% add additive gaussian white noise
dataSNR = 69; % in dB, set a stupidly high SNR to make them okay
space3Dnoised = awgn(space3D,dataSNR);
sg = reshape(space3Dnoised(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3Dnoised(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRnoised = snr2d(sg,bg);

% change the original value
space3Dori = space3D;
space3D = space3Dnoised;

% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
tic
space3Dbpf = filtfilt(b,a,space3D);
tbpf = toc;

organized_show_image(space3Dbpf)

%% denoise using SVD
tic;[denoisedSVDgamma0]  = svd_denoising(space3Dbpf,0); tSVDgamma0 = toc;
tic;[denoisedSVDgammahalf]  = svd_denoising(space3Dbpf,0.5); tSVDgammahalf = toc;
tic;[denoisedSVDgamma1]  = svd_denoising(space3Dbpf,1); tSVDgamma1 = toc;
tic;[denoisedPaper1EMDMI] = paper_denoising(space3Dbpf,'paper-1-emd-mi'); tPaper1EMDMI = toc;
tic;[denoisedPaper1Wavelet] = paper_denoising(space3Dbpf,'paper-1-wavelet'); tPaper1Wavelet = toc;
tic;[denoisedPaper2DWT] = paper_denoising(space3Dbpf,'paper-2-dwt'); tPaper2DWT = toc;
tic;[denoisedPaper2MODWT] = paper_denoising(space3Dbpf,'paper-2-modwt'); tPaper2MODWT = toc;

%% get evaluation metrics
evalBPF = get_metrics(space3Dbpf,space3Dori,tbpf);
evalSVDgamma0 = get_metrics(denoisedSVDgamma0,space3Dori,tSVDgamma0);
evalSVDgammahalf = get_metrics(denoisedSVDgammahalf,space3Dori,tSVDgammahalf);
evalSVDgamma1 = get_metrics(denoisedSVDgamma1,space3Dori,tSVDgamma1);
evalPaper1EMDMI = get_metrics(denoisedPaper1EMDMI,space3Dori,tPaper1EMDMI);
evalPaper1Wavelet = get_metrics(denoisedPaper1Wavelet,space3Dori,tPaper1Wavelet);
evalPaper2DWT = get_metrics(denoisedPaper2DWT,space3Dori,tPaper2DWT);
evalPaper2MODWT = get_metrics(denoisedPaper2MODWT,space3Dori,tPaper2MODWT);

%% save
savename = sprintf("complete_comparison_flated_v3/rounded_eval_%d_%d.mat",dataSNR,idx);
%save(savename,"-regexp","^SNR","^eval","dataSNR")

%% calculate C-mode diameter
organized_show_image(space3Dbpf);
organized_show_image(denoisedSVDgammahalf);
organized_show_image(denoisedPaper1EMDMI);
organized_show_image(denoisedPaper1Wavelet);
%% visualization
% organized_show_image(space3D);
% organized_show_image(space3Dbpf);
% organized_show_image(denoisedSVDgammahalf);
% organized_show_image(denoisedPaper1Wavelet);
% organized_show_image(denoisedPaper2DWT);
% organized_show_image(denoisedPaper2MODWT);
% organized_show_image(denoisedPaper1EMDMI);
% 
% 
% % one signal show
% figure
% subplot(241);plot(space3Dori(:,102,33));title("GT");xlim([0 1024])
% subplot(242);plot(space3D(:,102,33));title("noisy");xlim([0 1024])
% subplot(243);plot(space3Dbpf(:,102,33));title("bpf");xlim([0 1024])
% subplot(244);plot(denoisedSVDgamma1(:,102,33));title("svd");xlim([0 1024])
% subplot(245);plot(denoisedPaper1Wavelet(:,102,33));title("sym6-sure");xlim([0 1024])
% subplot(246);plot(denoisedPaper2DWT(:,102,33));title("sym4-univ.ths.");xlim([0 1024])
% subplot(247);plot(denoisedPaper2MODWT(:,102,33));title("sym4-modwt-univ.ths.");xlim([0 1024])
% subplot(248);plot(denoisedPaper1EMDMI(:,102,33));title("emd-mi");xlim([0 1024])

%% evaluation function
function [a] = get_metrics(D,ref,timerecord)
    a.cnr = get_cnr(D);
    a.snr = get_snr(D);
    a.ssim_cmode = get_ssim_cmode(D,ref);
    a.psnr_cmode = get_psnr_cmode(D,ref);
    a.time = timerecord;
end

function [a] = get_cnr(D)
    global_param;
    sgx = sgx_flat2; 
    sgy = sgy_flat2;
    bgx = bgx_flat2; 
    bgy = bgy_flat2;

    hb = abs(hilbert(D));
    cmode = squeeze(max(hb));
    sgmu = mean(cmode(sgx,sgy),"all");
    bgmu = mean(cmode(bgx,bgy),"all");
    sgsd = std(cmode(sgx,sgy),0,"all");
    bgsd = std(cmode(bgx,bgy),0,"all");
    a = (sgmu - bgmu) / (sqrt((sgsd^2 + bgsd^2) / 2));
    a = mag2db(a);
end


function [a] = get_snr(D)
    global_param;
    sgx = sgx_flat2; 
    sgy = sgy_flat2;
    bgx = bgx_flat2; 
    bgy = bgy_flat2;

    sg = reshape(D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
    bg = reshape(D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));

    a = snr2d(sg,bg);
end

function [a] = get_ssim_cmode(D,ref)
    % D
    hb = abs(hilbert(D));
    cmodeD = squeeze(max(hb));
    cmodeD = (cmodeD - min(cmodeD,[],"all")) / range(cmodeD,'all');

    % ref
    hb = abs(hilbert(ref));
    cmodeRef = squeeze(max(hb));
    cmodeRef = (cmodeRef - min(cmodeRef,[],"all")) / range(cmodeRef,"all");

    a = ssim(cmodeD,cmodeRef);
end

function [a] = get_psnr_cmode(D,ref)
    % D
    hb = abs(hilbert(D));
    cmodeD = squeeze(max(hb));
    cmodeD = (cmodeD - min(cmodeD,[],"all")) / range(cmodeD,'all');

    % ref
    hb = abs(hilbert(ref));
    cmodeRef = squeeze(max(hb));
    cmodeRef = (cmodeRef - min(cmodeRef,[],"all")) / range(cmodeRef,"all");

    a = psnr(cmodeD,cmodeRef);
end

function [outp] = snr2d(a,b)
    sz = size(a,2);
    outp = 0;
    for i =1:sz
        outp = outp + snr(a(:,i),b(:,i));
    end
    outp = outp / sz;
end