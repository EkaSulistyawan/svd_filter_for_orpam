close all
clear all

load("../../SharingPoint/data/cellular/compiled/17122022_USAF.mat")
space3D = space3D(:,1:200,1:200);
Fs = 5e9;

% Compute base SNR
global_param;
sgx = sgx_flat;
sgy = sgy_flat;
bgx = bgx_flat;
bgy = bgy_flat;

sg = reshape(space3D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRinit = snr(sg,bg);

% add additive gaussian white noise
dataSNR = 60; % in dB, set a stupidly high SNR to make them okay
space3Dnoised = awgn(space3D,dataSNR);
sg = reshape(space3Dnoised(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3Dnoised(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRnoised = snr(sg,bg);

% change the original value
space3Dori = space3D;
space3D = space3Dnoised;

% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3Dnoised);

%organized_show_image(space3Dbpf)

%% denoise using SVD
tic;[denoisedSVD]  = svd_denoising(space3Dbpf,0.5);toc;
tic;[denoisedPaper1Wavelet] = waveletfam_denoising(space3Dbpf,'paper-1-wavelet');toc;
tic;[denoisedPaper2DWT] = waveletfam_denoising(space3Dbpf,'paper-2-dwt');toc
tic;[denoisedPaper2MODWT] = waveletfam_denoising(space3Dbpf,'paper-2-modwt');toc

%%
organized_show_image(space3D);
organized_show_image(space3Dbpf);
organized_show_image(denoisedSVD);
organized_show_image(denoisedPaper1Wavelet);
organized_show_image(denoisedPaper2DWT);
organized_show_image(denoisedPaper2MODWT);
%% Do EMD on one of their signal
% signal = space3D(:,101,32);
% emd(signal);
% [imfs,residual] = emd(signal);
% figure;plot(imfs)

