clear all
close all

%% load data
load("../../SharingPoint/data/cellular/compiled/ori_100avg_4d_sphere.mat");

%%
space4D = space4D(:,:,1:100,1:100);

%%
avgnum = 100;
space3Dref = squeeze(mean(space4D(:,1:avgnum,:,:),2));

%%
avgnum = 1;
space3D = squeeze(mean(space4D(:,1:avgnum,:,:),2));

%%
% for nonflat
global_param;
sgx = sgx_varavg; 
sgy = sgy_varavg;
bgx = bgx_varavg; 
bgy = bgy_varavg;

sg = reshape(space3D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRinit = snr2d(sg,bg);

% BPF
Fs = 5e9;
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
tic
space3Dbpf = filtfilt(b,a,space3D);
tbpf = toc;

% denoise
tic;[denoisedSVDgamma0]  = svd_denoising(space3Dbpf,0); tSVDgamma0 = toc;
tic;[denoisedSVDgammahalf]  = svd_denoising(space3Dbpf,0.5); tSVDgammahalf = toc;
tic;[denoisedSVDgamma1]  = svd_denoising(space3Dbpf,1); tSVDgamma1 = toc;
tic;[denoisedPaper1EMDMI] = paper_denoising(space3Dbpf,'paper-1-emd-mi'); tPaper1EMDMI = toc;
tic;[denoisedPaper1Wavelet] = paper_denoising(space3Dbpf,'paper-1-wavelet'); tPaper1Wavelet = toc;
tic;[denoisedPaper2DWT] = paper_denoising(space3Dbpf,'paper-2-dwt'); tPaper2DWT = toc;
tic;[denoisedPaper2MODWT] = paper_denoising(space3Dbpf,'paper-2-modwt'); tPaper2MODWT = toc;

tic;[denoisedSVDgamma0_wobpf]  = svd_denoising(space3D,0); tSVDgamma0_wobpf = toc;
tic;[denoisedSVDgammahalf_wobpf]  = svd_denoising(space3D,0.5); tSVDgammahalf_wobpf = toc;
tic;[denoisedSVDgamma1_wobpf]  = svd_denoising(space3D,1); tSVDgamma1_wobpf = toc;
tic;[denoisedPaper1EMDMI_wobpf] = paper_denoising(space3D,'paper-1-emd-mi'); tPaper1EMDMI_wobpf = toc;
tic;[denoisedPaper1Wavelet_wobpf] = paper_denoising(space3D,'paper-1-wavelet'); tPaper1Wavelet_wobpf = toc;
tic;[denoisedPaper2DWT_wobpf] = paper_denoising(space3D,'paper-2-dwt'); tPaper2DWT_wobpf = toc;
tic;[denoisedPaper2MODWT_wobpf] = paper_denoising(space3D,'paper-2-modwt'); tPaper2MODWT_wobpf = toc;

%% evaluate

evalBPF = get_metrics(space3Dbpf,space3Dref,tbpf);
evalSVDgamma0 = get_metrics(denoisedSVDgamma0,space3Dref,tSVDgamma0);
evalSVDgammahalf = get_metrics(denoisedSVDgammahalf,space3Dref,tSVDgammahalf);
evalSVDgamma1 = get_metrics(denoisedSVDgamma1,space3Dref,tSVDgamma1);
evalPaper1EMDMI = get_metrics(denoisedPaper1EMDMI,space3Dref,tPaper1EMDMI);
evalPaper1Wavelet = get_metrics(denoisedPaper1Wavelet,space3Dref,tPaper1Wavelet);
evalPaper2DWT = get_metrics(denoisedPaper2DWT,space3Dref,tPaper2DWT);
evalPaper2MODWT = get_metrics(denoisedPaper2MODWT,space3Dref,tPaper2MODWT);

evalBPF_wobpf = get_metrics(space3Dbpf,space3Dref,tbpf);
evalSVDgamma0_wobpf = get_metrics(denoisedSVDgamma0_wobpf,space3Dref,tSVDgamma0_wobpf);
evalSVDgammahalf_wobpf = get_metrics(denoisedSVDgammahalf_wobpf,space3Dref,tSVDgammahalf_wobpf);
evalSVDgamma1_wobpf = get_metrics(denoisedSVDgamma1_wobpf,space3Dref,tSVDgamma1_wobpf);
evalPaper1EMDMI_wobpf = get_metrics(denoisedPaper1EMDMI_wobpf,space3Dref,tPaper1EMDMI_wobpf);
evalPaper1Wavelet_wobpf = get_metrics(denoisedPaper1Wavelet_wobpf,space3Dref,tPaper1Wavelet_wobpf);
evalPaper2DWT_wobpf = get_metrics(denoisedPaper2DWT_wobpf,space3Dref,tPaper2DWT_wobpf);
evalPaper2MODWT_wobpf = get_metrics(denoisedPaper2MODWT_wobpf,space3Dref,tPaper2MODWT_wobpf);

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
    sgx = sgx_round; 
    sgy = sgy_round;
    bgx = bgx_round; 
    bgy = bgy_round;

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
    sgx = sgx_round; 
    sgy = sgy_round;
    bgx = bgx_round; 
    bgy = bgy_round;

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