%% add new data for artifact removal
clear all
close all
load("../../SharingPoint/data/cellular/compiled/23022023_rbc.mat")
% load("../../Tohoku/Research Project/CS/Matlab/data/cellular/compiled/23022023_rbc.mat")
% show image
sz = size(space3D,2)-1;
space3D = space3D(:,1:sz,1:sz);
cmodefunc = @(x)(squeeze(max(abs(hilbert(x)))));
space3D2 = space3Draw;
%% modify
modx1 = 7:10   ; mody1 = 78:80; refx1 = 7; refy1 = 78;
modx2 = 151:156; mody2 = 48:52; refx2 = 151; refy2 = 48;

% dat2d1 = reshape(space3D(:,modx1,mody1),1024,numel(modx1)*numel(mody1));
% min1 = min(space3D(:,refx1,refy1));
% ran1 = range(space3D(:,refx1,refy1));

for i=1:numel(modx1)
    for j=1:numel(mody1)
        refsglvl = range(space3D(:,refx1,refy1));
        tgtsg = space3D2(:,modx1(i),mody1(j));
        tgtsglvl = range(tgtsg);

        space3D2(:,modx1(i),mody1(j)) = space3D2(:,modx1(i),mody1(j))*refsglvl/tgtsglvl;
    end
end

for i=1:numel(modx2)
    for j=1:numel(mody2)
        refsglvl = range(space3D(:,refx2,refy2));
        tgtsg = space3D2(:,modx2(i),mody2(j));
        tgtsglvl = range(tgtsg);

        space3D2(:,modx2(i),mody2(j)) = space3D2(:,modx2(i),mody2(j))*refsglvl/tgtsglvl;
    end
end
space3D = space3D2; 
im = cmodefunc(space3D);
figure;imagesc(im);axis image;colormap hot
%%
Fs = 5e9;

% for nonflat
global_param;
sgx = sgx_rbcArt; 
sgy = sgy_rbcArt;
bgx = bgx_rbcArt; 
bgy = bgy_rbcArt;

% sg = reshape(space3D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
% bg = reshape(space3D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
% SNRinit = snr2d(sg,bg);


get_snr(space3D)
get_cnr(space3D)


% BPF
fcl = 20e6;
fch = 80e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
tic
space3Dbpf = filtfilt(b,a,space3D);
tbpf = toc;
%%
D = space3D;
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[u,s,v] = svd(dat2d,"econ");

D = space3Dbpf;
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[ubpf,sbpf,vbpf] = svd(dat2d,"econ");

% save("rbc_artifact_USV.mat")
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

% get evaluation metrics
% space3Dori = space3D;
evalBPF = get_metrics(space3Dbpf,space3D,tbpf);
% evalSVDgamma0 = get_metrics(denoisedSVDgamma0,space3Dbpf,tSVDgamma0);
% evalSVDgammahalf = get_metrics(denoisedSVDgammahalf,space3Dbpf,tSVDgammahalf);
% evalSVDgamma1 = get_metrics(denoisedSVDgamma1,space3Dbpf,tSVDgamma1);
% evalPaper1EMDMI = get_metrics(denoisedPaper1EMDMI,space3Dbpf,tPaper1EMDMI);
evalPaper1Wavelet = get_metrics(denoisedPaper1Wavelet,space3D,tPaper1Wavelet);
evalPaper2DWT = get_metrics(denoisedPaper2DWT,space3D,tPaper2DWT);
evalPaper2MODWT = get_metrics(denoisedPaper2MODWT,space3D,tPaper2MODWT);
evalSVDwV = get_metrics(denoisedSVDwV,space3D,tSVDwV);
evalSVDwU = get_metrics(denoisedSVDwU,space3D,tSVDwU);
evalSVDwVnbpf = get_metrics(denoisedSVDwVnbpf,space3D,tSVDwV);
evalSVDwUnbpf = get_metrics(denoisedSVDwUnbpf,space3D,tSVDwU);

%
save("rbc_artifact3.mat")

%%

% im = cmodefunc(denoisedSVDgamma0);
% figure;imagesc(im);axis image;colormap hot;title("\gamma=0");axis off
% im = cmodefunc(denoisedSVDgammahalf);
% figure;imagesc(im);axis image;colormap hot;title("\gamma=0.5");axis off
% im = cmodefunc(denoisedSVDgamma1);
% figure;imagesc(im);axis image;colormap hot;title("\gamma=1");axis off
% im = cmodefunc(denoisedPaper1EMDMI);
% figure;imagesc(im);axis image;colormap hot;title("EMDMI");axis off
% im = cmodefunc(denoisedPaper1Wavelet);
% figure;imagesc(im);axis image;colormap hot;title("Wavelet");axis off
% im = cmodefunc(denoisedPaper2DWT);
% figure;imagesc(im);axis image;colormap hot;title("DWT");axis off
% im = cmodefunc(denoisedPaper2MODWT);
% figure;imagesc(im);axis image;colormap hot;title("MODWT");axis off

%% residual
% save("saved_others\rbc_artifact_test.mat")
% im = cmodefunc(space3D2 - denoisedSVDgamma0);
% figure;imagesc(im);axis image;colormap hot;title("\gamma=0");axis off
% im = cmodefunc(space3D2 - denoisedSVDgammahalf);
% figure;imagesc(im);axis image;colormap hot;title("\gamma=0.5");axis off
% im = cmodefunc(space3D2 - denoisedSVDgamma1);
% figure;imagesc(im);axis image;colormap hot;title("\gamma=1");axis off 
% im = cmodefunc(space3D2 - denoisedPaper1EMDMI);
% figure;imagesc(im);axis image;colormap hot;title("EMDMI");axis off
% im = cmodefunc(space3D2 - denoisedPaper1Wavelet);
% figure;imagesc(im);axis image;colormap hot;title("Wavelet");axis off
% im = cmodefunc(space3D2 - denoisedPaper2DWT);
% figure;imagesc(im);axis image;colormap hot;title("DWT");axis off
% im = cmodefunc(space3D2 - denoisedPaper2MODWT);
% figure;imagesc(im);axis image;colormap hot;title("MODWT");axis off

%% evaluation function

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
    sgx = sgx_rbcArt; 
    sgy = sgy_rbcArt;
    bgx = bgx_rbcArt; 
    bgy = bgy_rbcArt;

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
    sgx = sgx_rbcArt; 
    sgy = sgy_rbcArt;
    bgx = bgx_rbcArt; 
    bgy = bgy_rbcArt;

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

function [a] = range(p)
    a = max(p,[],"all") - min(p,[],'all');
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