%% load data 
load("../../SharingPoint/data/cellular/compiled/17122022_USAF.mat")
space3D = space3D(:,1:200,1:200);
Fs = 5e9;

% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);

% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,200*200);
[u,sflat,v] = svd(dat2d,"econ");

%% load round
load("../../SharingPoint/data/cellular/compiled/09102022_sphere.mat")
%load("../data/cellular/compiled/29112022_10umsphere.mat"); % same dimension

space3D = space3D(:,1:200,1:200);
Fs = 5e9;

% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);
hb = abs(hilbert(space3Dbpf));
cmode = squeeze(max(hb));
figure;imagesc(cmode);axis image;colormap gray

% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,200*200);
[u,sround,v] = svd(dat2d,"econ");

%%
figure;
yyaxis left
plot(diag(sflat));
yyaxis right
plot(diag(sround));
xlim([0 100])