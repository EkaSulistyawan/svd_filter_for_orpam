function [rec3d] = waveletfam_denoising(D,paperNum)
% THIS ARE COLLECTION OF PAPER IMPLEMENTED TECHNIQUE FOR DENOISING
% re-arrange
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);

% % 
% Let's try to replicate what have been reported by this paper:
% A Noise Reduction Method for Photoacoustic Imaging In Vivo Based on EMD 
% and Conditional Mutual Information
% 
% I was still skeptical on how they implement this technique properly,
% there was no proper code being shared nor exact equation.

if(paperNum == "paper-1-wavelet")
    disp("Run PAPER 1 WAVELET")
    [denoised] = wdenoise(dat2d,wmaxlev(1024,'sym6'),'Wavelet','sym6','DenoisingMethod','SURE');


% 
% Paper:
% Automated wavelet denoising of photoacoustic signals for circulating
% melanoma cell detection and burn image reconstruction
%
% There was no apparent reason in using MODWT since our data is 2^j
% so it's actually the same with a referenced implementation in their paper
% 
% the paper is:
% In vivo port wine stain depth determination using a photoacoustic probe
% 
% 

elseif(paperNum == "paper-2-dwt")
    disp("Run PAPER 2 DWT")
    [denoised] = wdenoise(dat2d,wmaxlev(1024,'sym4'),'Wavelet','sym4','DenoisingMethod','UniversalThreshold');
elseif(paperNum == "paper-2-modwt")
    disp("Run PAPER 2 MODWT")
    denoised = zeros(size(dat2d));
    for i=1:size(dat2d,2)
        signal = dat2d(:,i);
        [outp] = wden(signal,'modwtsqtwolog','s','mln',4,'sym4');
        denoised(:,i) = outp;
    end
end
rec3d= reshape(denoised,1024,imsize,imsize);
end

