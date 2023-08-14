%%
% Evaluation with Synthetic Noises
% Load the flat-surfaced data

close all
clear all

load("../../SharingPoint/data/cellular/compiled/17122022_USAF.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:200,1:200);
Fs = 5e9;

% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);


%% do filtering
[denoised]  = svd_denoising(space3Dbpf,0.5);

%% show image
organized_show_image(space3Dbpf)
organized_show_image(denoised)

%%                                                                          METRICS!!!
%% CNR
global_param;
sgx = sgx_flat;
sgy = sgy_flat;
bgx = bgx_flat;
bgy = bgy_flat;

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