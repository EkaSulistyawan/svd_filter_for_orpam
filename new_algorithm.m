% %
% we rely on the energy of the signal
% in the sense that PA signal is a wavelet --> 
% occupy finite energy in both frequency and time domain
% 
% READ MALLAT's BOOK!!!!
% %

%% load data
load("../data/cellular/compiled/13122022_25uWRBC.mat");
imsize = size(space3D,2) - mod(size(space3D,2),10);

space3D = space3D(:,1:imsize,1:imsize);
Fs = 5e9;

% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);

% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,imsize^2);
[u,s,v] = svd(dat2d,"econ");

%% energy based U 
sig = u;
minsig = repmat(min(sig),[size(u,2) 1]);
rangesig = repmat(range(sig),[size(u,2) 1]);
sig = (sig - minsig) ./ rangesig;
sig = (sig - 0.5)*2;
energy = -sum(sig.^2);
energy = (energy - min(energy))./ range(energy);

% multiplier
w0 = energy' .* diag(s);
%wnotselected = (1-energy)' .* diag(s);

%% Objectness based V
v3d = reshape(v,imsize,imsize,1024);
objectness = get_dctweight(v3d);


% multiplier
w1 = objectness' .* diag(s);
%wnotselected = (1-energy)' .* diag(s);

%% reconstruct 
gamma = 1; % range 0 ~ 1
rec = u*diag((gamma*w0) + (1-gamma)*w1)*v';
rec3d= reshape(rec,1024,imsize,imsize);

% rem = u*diag(wnotselected)*v';
% rem3d= reshape(rem,1024,imsize,imsize);

%% show image
hb = abs(hilbert(space3Dbpf));
cmode = squeeze(max(hb));
minCaxis = min(cmode,[],"all");
maxCaxis = max(cmode,[],"all");

figure;imagesc(cmode);axis image;colormap hot;%caxis([minCaxis maxCaxis])
figure;imagesc(squeeze(hb(:,:,1)));colormap hot

hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
figure;imagesc(cmode);axis image;colormap hot;%caxis([minCaxis maxCaxis])
figure;imagesc(squeeze(hb(:,:,1)));colormap hot

% hb = abs(hilbert(rem3d));
% cmode = squeeze(max(hb));
% figure;imagesc(cmode);axis image;colormap hot;%caxis([minCaxis maxCaxis])
% figure;imagesc(squeeze(hb(:,:,1)));colormap hot