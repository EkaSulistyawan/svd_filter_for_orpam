% Photoacoustic Signal Extraction by Singular Value Decomposition with Weighting Matrix
% Author: I Gede Eka Sulistyawan (sulis347@gmail.com), 2023
% 
% Usage
% rec3d = svd_denoising_rank(D,gamma)
% 
% D    : input 3D signal (Nt, Nx, Ny), with A-lines at 1st dimension
% gamma: continuous value [0~1], control over strength of weight (wU or wV)
%           gamma=0  , only utilize wV
%           gamma=0.5, utilize half wU and half wV
%           gamma=1  , only utilize wU
% rec3D: the 3D cleaned signal

% referring to the paper
% P = UP SP VP'
% P' = UP (SP*(wU or wV)) VP'
% 
% P is reshaped D
% rec3D is reshaped P'
% gamma controls over wU or wV
% 

function [rec3d,wU,wVnorm,wVraw] = svd_denoising_ius(D,gamma)

% re-arrange
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[u,s,v] = svd(dat2d,"econ");

% wU
uu = normalize(u,1,"range")-0.5;
uu = sum(uu.^2);
wU = normalize(-uu,'range');
wUraw = wU;
wU = wU' .* diag(s);

% wV
v3d = reshape(v,imsize,imsize,1024);
ss = pagesvd(v3d,'econ');
vv = sum(squeeze(ss));
wV = normalize(-vv,'range');
wVraw = vv;
wVnorm = wV;
wV = wV' .* diag(s);


%rec = u*diag(wU.*wV)*v';
rec = u*diag(wV)*v';
rec3d= reshape(rec,1024,imsize,imsize);


end

