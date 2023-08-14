%%                                                                          SIMULATED NOISES
% This technique is dependant on SNRraw, so make sure you compute them
% first. We will use awgn to simulate the noise level, i.e., decrease the
% SNR appropriately
close all
clear all

load("../../SharingPoint/data/cellular/compiled/17122022_USAF.mat")
%load("../../SharingPoint/data/cellular/compiled/09102022_sphere.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
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
addNoise = 55; % in dB
space3Dnoised = awgn(space3D,addNoise);
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


%% do filtering
[denoised]  = svd_denoising(space3Dbpf,0.5);

%% show image
organized_show_image(space3Dori)
organized_show_image(space3D)
organized_show_image(space3Dbpf)
organized_show_image(denoised)

%%                                                                          METRICS!!!
%% CNR

hb = abs(hilbert(space3D));
cmode = squeeze(max(hb));
sgmu = mean(cmode(sgx,sgy),"all");
bgmu = mean(cmode(bgx,bgy),"all");
sgsd = std(cmode(sgx,sgy),0,"all");
bgsd = std(cmode(bgx,bgy),0,"all");
CNRraw = (sgmu - bgmu) / (sqrt((sgsd^2 + bgsd^2) / 2));
CNRraw = mag2db(CNRraw)

hb = abs(hilbert(space3Dbpf));
cmode = squeeze(max(hb));
sgmu = mean(cmode(sgx,sgy),"all");
bgmu = mean(cmode(bgx,bgy),"all");
sgsd = std(cmode(sgx,sgy),0,"all");
bgsd = std(cmode(bgx,bgy),0,"all");
CNRbpf = (sgmu - bgmu) / (sqrt((sgsd^2 + bgsd^2) / 2));
CNRbpf = mag2db(CNRbpf)

hb = abs(hilbert(denoised));
cmode = squeeze(max(hb));
sgmu = mean(cmode(sgx,sgy),"all");
bgmu = mean(cmode(bgx,bgy),"all");
sgsd = std(cmode(sgx,sgy),0,"all");
bgsd = std(cmode(bgx,bgy),0,"all");
CNRsvd = (sgmu - bgmu) / (sqrt((sgsd^2 + bgsd^2) / 2));
CNRsvd = mag2db(CNRsvd)

%% SNR
sg = reshape(space3D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRraw = snr(sg,bg)

sg = reshape(space3Dbpf(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3Dbpf(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRbpf = snr(sg,bg)

sg = reshape(denoised(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(denoised(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRsvd = snr(sg,bg)

%% SSIM against ori
cmodefunc = @(x) squeeze(max(abs(hilbert(x))));
SSIMraw = ssim(cmodefunc(space3D),cmodefunc(space3Dori))
SSIMbpf = ssim(cmodefunc(space3Dbpf),cmodefunc(space3Dori))
SSIMsvd = ssim(cmodefunc(denoised),cmodefunc(space3Dori))
%%
PSNRraw = psnr(abs(hilbert(space3D)),abs(hilbert(space3Dori)))
PSNRbpf = psnr(abs(hilbert(space3Dbpf)),abs(hilbert(space3Dori)))
PSNRsvd = psnr(abs(hilbert(denoised)),abs(hilbert(space3Dori)))