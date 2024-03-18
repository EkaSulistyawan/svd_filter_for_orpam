%% Comparison only on one signal
clear all
close all


load("../../SharingPoint/data/cellular/compiled/07062023_rbcFull.mat")
%load("../../Tohoku/Research Project/CS/Matlab/data/cellular/compiled/07062023_rbcFull.mat")
%%
sz = size(space3D,2)-1;
space3D = space3D(:,1:sz,1:sz);
Fs = 5e9;

% for nonflat
global_param;
sgx = sgx_rbc; 
sgy = sgy_rbc;
bgx = bgx_rbc; 
bgy = bgy_rbc;

sg = reshape(space3D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));

get_snr(space3D)
get_cnr(space3D)

% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl*2/Fs fch*2/Fs]);
tic
space3Dbpf = filtfilt(b,a,space3D);
tbpf = toc;

cutoff = 1e6;
[b,a] = butter(4,cutoff / (Fs/2),'high');
space3Dhpf = filtfilt(b,a,space3D);

organized_show_image(space3D)
%% show 
D = space3D;
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[u,s,v] = svd(dat2d,"econ");

D = space3Dbpf;
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[ubpf,sbpf,vbpf] = svd(dat2d,"econ");

%save("rbc_denoising_USV.mat")
%%
% close all
% cmode_func = @(x)(squeeze(max(abs(hilbert(x)))));
% cmode = cmode_func(denoisedSVDgammahalf);
% figure;imagesc(cmode');axis image;colormap hot
% figure;plot(cmode(1:33,71)); %% for biconcave
% figure;plot(cmode(58:95,112)); %% for biconcave
% figure;plot(cmode(95:120,112)); %% for biconcave
% vfigure;plot(cmode(36:66,74)); %% for biconcave

%% denoise using SVD
% tic;[denoisedSVDgamma0]  = svd_denoising(space3Dbpf,0); tSVDgamma0 = toc;
% tic;[denoisedSVDgammahalf]  = svd_denoising(space3Dbpf,0.5); tSVDgammahalf = toc;
% tic;[denoisedSVDgamma1]  = svd_denoising(space3Dbpf,1); tSVDgamma1 = toc;
% tic;[denoisedPaper1EMDMI] = paper_denoising(space3Dbpf,'paper-1-emd-mi'); tPaper1EMDMI = toc;
tic;[denoisedPaper1Wavelet] = paper_denoising(space3D,'paper-1-wavelet'); tPaper1Wavelet = toc;
tic;[denoisedPaper2DWT] = paper_denoising(space3D,'paper-2-dwt'); tPaper2DWT = toc;
tic;[denoisedPaper2MODWT] = paper_denoising(space3D,'paper-2-modwt'); tPaper2MODWT = toc;
tic;[denoisedSVDwV]  = svd_denoising_rank(space3Dbpf,0); tSVDwV = toc;% wV
tic;[denoisedSVDwU]  = svd_denoising_rank(space3Dbpf,1); tSVDwU = toc; 

tic;[denoisedSVDwVnbpf]  = svd_denoising_rank(space3D,0); tSVDwV = toc;% wV
tic;[denoisedSVDwUnbpf]  = svd_denoising_rank(space3D,1); tSVDwU = toc;  

tic;[denoisedSVDwVlet]  = svd_denoising_rank(denoisedPaper2MODWT,0); tSVDwV = toc;% wV
tic;[denoisedSVDwUlet]  = svd_denoising_rank(denoisedPaper2MODWT,1); tSVDwU = toc; 

%% get evaluation metrics
evalBPF = get_metrics(space3Dbpf,space3D,tbpf);
% evalSVDgamma0 = get_metrics(denoisedSVDgamma0,space3Dbpf,tSVDgamma0);
% evalSVDgammahalf = get_metrics(denoisedSVDgammahalf,space3Dori,tSVDgammahalf);
% evalSVDgamma1 = get_metrics(denoisedSVDgamma1,space3Dbpf,tSVDgamma1);
% evalPaper1EMDMI = get_metrics(denoisedPaper1EMDMI,space3Dori,tPaper1EMDMI);
evalPaper1Wavelet = get_metrics(denoisedPaper1Wavelet,space3D,tPaper1Wavelet);
evalPaper2DWT = get_metrics(denoisedPaper2DWT,space3D,tPaper2DWT);
evalPaper2MODWT = get_metrics(denoisedPaper2MODWT,space3D,tPaper2MODWT);
%
evalSVDwV = get_metrics(denoisedSVDwV,space3D,tSVDwV+tbpf);
evalSVDwU = get_metrics(denoisedSVDwU,space3D,tSVDwU+tbpf);

evalSVDwVnbpf = get_metrics(denoisedSVDwVnbpf,space3D,tSVDwV);
evalSVDwUnbpf = get_metrics(denoisedSVDwUnbpf,space3D,tSVDwU);

evalSVDwVlet = get_metrics(denoisedSVDwVlet,space3D,tSVDwV);
evalSVDwUlet = get_metrics(denoisedSVDwUlet,space3D,tSVDwU);

%%
% save("rbc_denoising_pluswavelet.mat")
%% save
%savename = sprintf("complete_comparison_rounded_v3/rounded_eval_%d_%d.mat",dataSNR,idx);
%save(savename,"-regexp","^SNR","^eval","dataSNR")
% get_metrics(denoisedPaper1Wavelet,space3Dbpf,tSVDwV);
% get_metrics(denoisedPaper1Wavelet,space3D,tPaper1Wavelet)
% organized_show_image(denoisedPaper1Wavelet)
% %%
% organized_show_image(denoisedSVDgammahalf);
% organized_show_image(denoisedSVDgamma0);
% organized_show_image(denoisedPaper1EMDMI);
% organized_show_image(denoisedPaper1Wavelet);
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
function [a] = range(p)
    a = max(p,[],"all") - min(p,[],'all');
end

function [a,sd] = pagessim(D,ref)
    D = abs(hilbert(D));
    R = abs(hilbert(ref));

    kk = zeros(1,size(D,1));
    for i=1:size(D,1)
        Dslice = squeeze(D(i,:,:));
        Dslice = rescale(Dslice);

        Rslice = squeeze(R(i,:,:));
        Rslice = rescale(Rslice);
        kk(i) = ssim(Dslice,Rslice);
    end
    figure;plot(kk)
    a = mean(kk);
    sd = std(kk);
end

function [a] = get_metrics(D,ref,timerecord)
    a.snr = get_snr(D);
    a.cnr = get_cnr(D);
    a.time = timerecord;
    
    % [b,c] = pagessim(D,ref);
    % a.ssim_page = b;
    % a.ssim_page_sd = c;
    a.ssim_cmode = get_ssim_cmode(D,ref);
    a.psnr_cmode = get_psnr_cmode(D,ref);
    
end

function [a] = get_cnr(D)
    global_param;
    sgx = sgx_rbc; 
    sgy = sgy_rbc;
    bgx = bgx_rbc; 
    bgy = bgy_rbc;

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
    sgx = sgx_rbc; 
    sgy = sgy_rbc;
    bgx = bgx_rbc; 
    bgy = bgy_rbc;

    sg = reshape(D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
    bg = reshape(D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));

    a = snr2d(sg,bg);
end

function [a] = get_ssim_cmode(D,ref)
    % D
    hb = abs(hilbert(D));
    cmodeD = squeeze(max(hb));
    cmodeD = (cmodeD - min(cmodeD,[],"all")) / range(cmodeD);

    % ref
    hb = abs(hilbert(ref));
    cmodeRef = squeeze(max(hb));
    cmodeRef = (cmodeRef - min(cmodeRef,[],"all")) / range(cmodeRef);

    a = ssim(cmodeD,cmodeRef);
end

function [a] = get_psnr_cmode(D,ref)
    % D
    hb = abs(hilbert(D));
    cmodeD = squeeze(max(hb));
    cmodeD = (cmodeD - min(cmodeD,[],"all")) / range(cmodeD);

    % ref
    hb = abs(hilbert(ref));
    cmodeRef = squeeze(max(hb));
    cmodeRef = (cmodeRef - min(cmodeRef,[],"all")) / range(cmodeRef);

    a = psnr(cmodeD,cmodeRef);
end

function [outp] = snr2d(a,b)
    sz = size(a,2);
    outp = 0;
    for i =1:sz
        outp = outp + custom_snr(a(:,i),b(:,i));
    end
    outp = outp / sz;
end

function [outp] = custom_snr(a,b)
    lowlim = 10e6;
    hghlim = 100e6;

    global_param

    faxis = (1:512)*5e9/1024;

    idx = find((faxis > lowlim) & (faxis < hghlim));
    fta = abs(fft(a));
    ftb = abs(fft(b));
    maxdbsg = max(fta(idx));
    maxdbns = max(ftb(idx));

    outp = mag2db(maxdbsg/maxdbns);
end
