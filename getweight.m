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

function [wU,wV,ss] = getweight(u,v,s)

% wU
uu = normalize(u,1,"range")-0.5;
uu = sum(uu.^2);
wU = normalize(-uu,'range');
wU = wU' ;

% wV
v3d = reshape(v,imsize,imsize,1024);
ss = pagesvd(v3d,'econ');
vv = sum(squeeze(ss));
wV = normalize(-vv,'range');
wV = wV' ;
ss = diag(s);

% rec = u*diag((gamma*wU) + (1-gamma)*wV)*v';
% rec3d= reshape(rec,1024,imsize,imsize);


end

