%% load data 
close all
clear all
%load("../data/cellular/compiled/11102022_sphere2.mat");
load("../../SharingPoint/data/cellular/compiled/09102022_sphere.mat")
imsize = size(space3D,2);
space3D = space3D(:,1:imsize,1:imsize);
Fs = 5e9;

%% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);

%% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,imsize*imsize);
[u,s,v] = svd(dat2d,"econ");

figure;imagesc(dat2d);axis image;colormap jet

%% Individual; get the first C-mode, B-mode and the signal
sel = [1:10];
% see the A line
subplot(3,1,1);plot(sum(u(:,sel),2));
xlim([1, 1024])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,imsize,imsize);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,1,2);imagesc(cmode);axis off;colormap hot

bmode = squeeze(hb(:,99,:));
subplot(3,1,3);imagesc(bmode);axis off;colormap hot



%% Collection; get the first C-mode, B-mode and the signal
numSel = 6;
axesCount = 1;

figure;
for p=1:numSel
sel = p;
% see the A line
maxval = max(u(:,1:numSel),[],"all");
minval = min(u(:,1:numSel),[],"all");
subplot(3,numSel,axesCount);plot(u(:,sel));
ylim([minval, maxval])
xlim([1, 1024])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,imsize,imsize);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,numSel,axesCount+numSel);imagesc(cmode);axis off;colormap hot

bmode = squeeze(hb(:,36,:));
subplot(3,numSel,axesCount+2*numSel);imagesc(bmode);axis off;colormap hot

axesCount = axesCount+1;
end

%% Collection; get the first C-mode, B-mode and the signal
numSel = 6;
axesCount = 1;

figure;
for p=1:numSel
sel = p + numSel;
% see the A line
maxval = max(u(:,1:numSel),[],"all");
minval = min(u(:,1:numSel),[],"all");
subplot(3,numSel,axesCount);plot(u(:,sel));
ylim([minval, maxval])
xlim([1, 1024])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,imsize,imsize);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,numSel,axesCount+numSel);imagesc(cmode);axis off;colormap hot

bmode = squeeze(hb(:,36,:));
subplot(3,numSel,axesCount+2*numSel);imagesc(bmode);axis off;colormap hot

axesCount = axesCount+1;
end

%% Collection; get the first C-mode, B-mode and the signal
numSel = 6;
axesCount = 1;

figure;
for p=1:numSel
sel = p + 2*numSel;
% see the A line
maxval = max(u(:,1:numSel),[],"all");
minval = min(u(:,1:numSel),[],"all");
subplot(3,numSel,axesCount);plot(u(:,sel));
ylim([minval, maxval])
xlim([1, 1024])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,imsize,imsize);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,numSel,axesCount+numSel);imagesc(cmode);axis off;colormap hot

bmode = squeeze(hb(:,36,:));
subplot(3,numSel,axesCount+2*numSel);imagesc(bmode);axis off;colormap hot

axesCount = axesCount+1;
end