
%% load data 
close all
clear all
load("../../SharingPoint/data/cellular/compiled/09102022_sphere.mat")
space3D = space3D(:,1:200,1:200);
Fs = 5e9;

% for nonflat
global_param;
sgx = sgx_round; 
sgy = sgy_round;
bgx = bgx_round; 
bgy = bgy_round;

sg = reshape(space3D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
bg = reshape(space3D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
SNRinit = snr(sg,bg);

% add additive gaussian white noise
dataSNR = 50; % in dB, set a stupidly high SNR to make them okay
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
tic
space3Dbpf = filtfilt(b,a,space3Dnoised);
tbpf = toc;
%%
imsize=200;
% SVD on bandpassed filter data
dat2d = reshape(space3D,1024,imsize*imsize);
[uraw,sraw,vraw] = svd(dat2d,"econ");

dat2d = reshape(space3Dbpf,1024,imsize*imsize);
[ubpf,sbpf,vbpf] = svd(dat2d,"econ");


%% with SVD
taxis = (1:1024)*1e9/Fs;
sg = space3D(:,100,150);
emdraw = emd(sg);

sg = space3Dbpf(:,100,150);
emdbpf = emd(sg);

%% 
close all
figure("Position",[3251,816,486,226]);
plot(taxis,space3D(:,100,150))
title("Raw Signal",'FontName','Times New Roman')
xlabel('Time (s)','FontName','Times New Roman','FontSize',12,'FontWeight','bold')
ylabel('Intensity (a.u.)','FontName','Times New Roman','FontSize',12,'FontWeight','bold')
set(gca,'FontName','Times New Roman','FontSize',12)
axis tight

ax = gcf;
exportgraphics(ax,'./exported_images/fig11_raw.emf')

figure("Position",[2045,108,1281,963]);
ax = tiledlayout(6,4,'TileSpacing','tight');
for i=1:6

nexttile()
psiv = zeros(size(sraw));
psiv(i,i) = 1;
dat3d = reshape(uraw*psiv*sraw*vraw',1024,200,200);
plot(taxis,dat3d(:,100,150))
ylabel(sprintf("Component \n%d",i))
if i==1;title("SVD without BPF",'FontName','Times New Roman');end
axis tight
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile()
psiv = zeros(size(sbpf));
psiv(i,i) = 1;
dat3d = reshape(ubpf*psiv*sbpf*vbpf',1024,200,200);
plot(taxis,dat3d(:,100,150))
if i==1;title("SVD with BPF",'FontName','Times New Roman');end
axis tight
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile()
plot(taxis,emdraw(:,i))
if i==1;title("EMD without BPF",'FontName','Times New Roman');end
axis tight
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile()
if(i < 4)
    plot(taxis,emdbpf(:,i))
    if i==1;title("EMD with BPF",'FontName','Times New Roman');end
    axis tight
else
    axis off
end
set(gca,'FontName','Times New Roman','FontSize',12)

end

xlabel(ax,'Time (s)','FontName','Times New Roman','FontSize',12,'FontWeight','bold')
ylabel(ax,'Intensity (a.u.)','FontName','Times New Roman','FontSize',12,'FontWeight','bold')

ax = gcf;
exportgraphics(ax,'./exported_images/fig11_compare.emf')