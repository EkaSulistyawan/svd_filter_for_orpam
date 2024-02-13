function [rec3d] = svd_denoising_rank(D,gamma)

% re-arrange
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[u,s,v] = svd(dat2d,"econ");

% wU
uu = normalize(u,1,"range")-0.5;
uu = sum(uu.^2);
wU = normalize(-uu,'range');
wU = wU' .* diag(s);

% wV
v3d = reshape(v,imsize,imsize,1024);
ss = pagesvd(v3d,'econ');
vv = sum(squeeze(ss));
wV = normalize(-vv,'range');
wV = wV' .* diag(s);


rec = u*diag((gamma*wU) + (1-gamma)*wV)*v';
rec3d= reshape(rec,1024,imsize,imsize);


end

