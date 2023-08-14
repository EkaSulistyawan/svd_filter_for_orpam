function [rec3d] = svd_denoising(D,gamma)

% re-arrange
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[u,s,v] = svd(dat2d,"econ");

sig = u;
minsig = repmat(min(sig),[size(u,2) 1]);
rangesig = repmat(range(sig),[size(u,2) 1]);
sig = (sig - minsig) ./ rangesig;
sig = (sig - 0.5)*2;
energy = -sum(sig.^2);
energy = (energy - min(energy))./ range(energy);

% multiplier
wU = energy' .* diag(s);

v3d = reshape(v,imsize,imsize,1024);
objectness = get_dctweight(v3d);


% multiplier
wV = objectness' .* diag(s);


rec = u*diag((gamma*wU) + (1-gamma)*wV)*v';
rec3d= reshape(rec,1024,imsize,imsize);


end

