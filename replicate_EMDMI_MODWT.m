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
dataSNR = 70; % in dB, set a stupidly high SNR to make them okay
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
D = space3Dbpf;
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
denoised = zeros(size(dat2d));
tic
parfor sgsel = 1:size(dat2d,2)
    signal = dat2d(:,sgsel);
    [imfs,residual] = emd(signal);
    
    partition = ceil(size(imfs,2)/2);
    tail = size(imfs,2);
    % find the partition by MI
    cond = true;
    while cond
        mi_with_noisy_signal = mi(signal,imfs(:,partition));
        selected_group = partition+1:tail;
        mi_with_selected_set = mi(imfs(:,partition),sum(imfs(:,selected_group),2));
        score = mi_with_noisy_signal - mi_with_selected_set;

        if(score < 0)% closer to group than noisy
            partition = partition -1;
        else
            cond = false;
        end

        if partition == 0
            cond = false;
        end
    end
    

    high_f_grp = 1:partition;
    low_f_grp = partition+1:tail;
    rec = sum(imfs(:,low_f_grp),2);
    % apply thresholding to high_f_grp
    for i =1:numel(high_f_grp)
        imf_t = imfs(:,high_f_grp(i));
        rec = rec + wthresh(imf_t,'h',thselect(imf_t,'rigrsure'));
    end
    % for the thresholding: https://www.mathworks.com/help/wavelet/ref/thselect.html

    denoised(:,sgsel) = rec;

end
toc
%%
rec3d= reshape(denoised,1024,imsize,imsize);
organized_show_image(rec3d)