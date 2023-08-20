%% compute avg 1
load("var_avg_data\splitted_indiv_avg_1.mat")
space3D = space3D(:,1:100,1:100);
Fs = 5e9;

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
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
tic
space3Dbpf = filtfilt(b,a,space3D);
tbpf = toc;


% denoise using SVD
tic;[denoisedSVDgamma0]  = svd_denoising(space3Dbpf,0); tSVDgamma0 = toc;
tic;[denoisedSVDgammahalf]  = svd_denoising(space3Dbpf,0.5); tSVDgammahalf = toc;
tic;[denoisedSVDgamma1]  = svd_denoising(space3Dbpf,1); tSVDgamma1 = toc;
tic;[denoisedPaper1EMDMI] = paper_denoising(space3Dbpf,'paper-1-emd-mi'); tPaper1EMDMI = toc;
tic;[denoisedPaper1Wavelet] = paper_denoising(space3Dbpf,'paper-1-wavelet'); tPaper1Wavelet = toc;
tic;[denoisedPaper2DWT] = paper_denoising(space3Dbpf,'paper-2-dwt'); tPaper2DWT = toc;
tic;[denoisedPaper2MODWT] = paper_denoising(space3Dbpf,'paper-2-modwt'); tPaper2MODWT = toc;

result_avg_1.raw = space3D;
result_avg_1.bpf = space3Dbpf;
result_avg_1.SVDgamma0 = denoisedSVDgamma0;
result_avg_1.SVDgammahalf = denoisedSVDgammahalf;
result_avg_1.SVDgamma1 = denoisedSVDgamma1;
result_avg_1.Paper1EMDMI = denoisedPaper1EMDMI;
result_avg_1.Paper1Wavelet = denoisedPaper1Wavelet;
result_avg_1.Paper2DWT = denoisedPaper2DWT;
result_avg_1.Paper2MODWT = denoisedPaper2MODWT;


%% compute avg 1
load("var_avg_data\splitted_indiv_avg_10.mat")
space3D = space3D(:,1:100,1:100);
Fs = 5e9;

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
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
tic
space3Dbpf = filtfilt(b,a,space3D);
tbpf = toc;


% denoise using SVD
tic;[denoisedSVDgamma0]  = svd_denoising(space3Dbpf,0); tSVDgamma0 = toc;
tic;[denoisedSVDgammahalf]  = svd_denoising(space3Dbpf,0.5); tSVDgammahalf = toc;
tic;[denoisedSVDgamma1]  = svd_denoising(space3Dbpf,1); tSVDgamma1 = toc;
tic;[denoisedPaper1EMDMI] = paper_denoising(space3Dbpf,'paper-1-emd-mi'); tPaper1EMDMI = toc;
tic;[denoisedPaper1Wavelet] = paper_denoising(space3Dbpf,'paper-1-wavelet'); tPaper1Wavelet = toc;
tic;[denoisedPaper2DWT] = paper_denoising(space3Dbpf,'paper-2-dwt'); tPaper2DWT = toc;
tic;[denoisedPaper2MODWT] = paper_denoising(space3Dbpf,'paper-2-modwt'); tPaper2MODWT = toc;

result_avg_10.raw = space3D;
result_avg_10.bpf = space3Dbpf;
result_avg_10.SVDgamma0 = denoisedSVDgamma0;
result_avg_10.SVDgammahalf = denoisedSVDgammahalf;
result_avg_10.SVDgamma1 = denoisedSVDgamma1;
result_avg_10.Paper1EMDMI = denoisedPaper1EMDMI;
result_avg_10.Paper1Wavelet = denoisedPaper1Wavelet;
result_avg_10.Paper2DWT = denoisedPaper2DWT;
result_avg_10.Paper2MODWT = denoisedPaper2MODWT;

%% compute avg 1
load("var_avg_data\splitted_indiv_avg_100.mat")
space3D = space3D(:,1:100,1:100);
Fs = 5e9;

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
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
tic
space3Dbpf = filtfilt(b,a,space3D);
tbpf = toc;


% denoise using SVD
tic;[denoisedSVDgamma0]  = svd_denoising(space3Dbpf,0); tSVDgamma0 = toc;
tic;[denoisedSVDgammahalf]  = svd_denoising(space3Dbpf,0.5); tSVDgammahalf = toc;
tic;[denoisedSVDgamma1]  = svd_denoising(space3Dbpf,1); tSVDgamma1 = toc;
tic;[denoisedPaper1EMDMI] = paper_denoising(space3Dbpf,'paper-1-emd-mi'); tPaper1EMDMI = toc;
tic;[denoisedPaper1Wavelet] = paper_denoising(space3Dbpf,'paper-1-wavelet'); tPaper1Wavelet = toc;
tic;[denoisedPaper2DWT] = paper_denoising(space3Dbpf,'paper-2-dwt'); tPaper2DWT = toc;
tic;[denoisedPaper2MODWT] = paper_denoising(space3Dbpf,'paper-2-modwt'); tPaper2MODWT = toc;

result_avg_100.raw = space3D;
result_avg_100.bpf = space3Dbpf;
result_avg_100.SVDgamma0 = denoisedSVDgamma0;
result_avg_100.SVDgammahalf = denoisedSVDgammahalf;
result_avg_100.SVDgamma1 = denoisedSVDgamma1;
result_avg_100.Paper1EMDMI = denoisedPaper1EMDMI;
result_avg_100.Paper1Wavelet = denoisedPaper1Wavelet;
result_avg_100.Paper2DWT = denoisedPaper2DWT;
result_avg_100.Paper2MODWT = denoisedPaper2MODWT;


clearvars -regexp ^denoised
save("var_avg_data\compiled_result.mat")


% %% get evaluation metrics
% evalBPF = get_metrics(space3Dbpf,space3Dori,tbpf);
% evalSVDgamma0 = get_metrics(denoisedSVDgamma0,space3Dori,tSVDgamma0);
% evalSVDgammahalf = get_metrics(denoisedSVDgammahalf,space3Dori,tSVDgammahalf);
% evalSVDgamma1 = get_metrics(denoisedSVDgamma1,space3Dori,tSVDgamma1);
% evalPaper1EMDMI = get_metrics(denoisedPaper1EMDMI,space3Dori,tPaper1EMDMI);
% evalPaper1Wavelet = get_metrics(denoisedPaper1Wavelet,space3Dori,tPaper1Wavelet);
% evalPaper2DWT = get_metrics(denoisedPaper2DWT,space3Dori,tPaper2DWT);
% evalPaper2MODWT = get_metrics(denoisedPaper2MODWT,space3Dori,tPaper2MODWT);





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
    sgx = sgx_flat; 
    sgy = sgy_flat;
    bgx = bgx_flat; 
    bgy = bgy_flat;

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
    sgx = sgx_flat; 
    sgy = sgy_flat;
    bgx = bgx_flat; 
    bgy = bgy_flat;

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