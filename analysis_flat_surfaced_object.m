%% load data 
load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:200,1:200);
Fs = 5e9;

%% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);

%% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,200*200);
[u,s,v] = svd(dat2d,"econ");

%% Individual; get the first C-mode, B-mode and the signal
sel = [1:10];
% see the A line
subplot(3,1,1);plot(sum(u(:,sel),2));
xlim([1, 1024])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,200,200);

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
rec3d= reshape(rec,1024,200,200);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,numSel,axesCount+numSel);imagesc(cmode);axis off;colormap hot

bmode = squeeze(hb(:,99,:));
subplot(3,numSel,axesCount+2*numSel);imagesc(bmode);axis off;colormap hot

axesCount = axesCount+1;
end