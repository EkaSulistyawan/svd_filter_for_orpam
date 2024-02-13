load("../../Tohoku/Research Project/CS/Matlab/data/cellular/compiled/09102022_sphere.mat")
sz = 200;
space3D = space3D(:,1:sz,1:sz);
Fs = 5e9;

space2D = reshape(space3D,1024,sz*sz);

% no filter

[up,sp,vp] = svd(space2D,'econ');

h = up(:,1);
T = 1024;
H = zeros(T,T);
for i=1:T
    H(i,:) = circshift(h,i);
end
H = circshift(H,512,2);

[uh,sh,vh] = svd(H);

% apply bpf
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
space2Df = filtfilt(b,a,space2D);
[upf,spf,vpf] = svd(space2Df,'econ');

h = upf(:,1);
T = 1024;
H = zeros(T,T);
for i=1:T
    H(i,:) = circshift(h,i);
end
H = circshift(H,512,2);

[uhf,shf,vhf] = svd(H);

%% Proof 1: some U are factorized signal over depth and its V correspondingly

%no filter 
figure;
nr = 6;nc=6;
tiledlayout(nr,nc,'TileSpacing','tight','Padding','none');

for ii=1:nr*nc
    nexttile
    imagesc(reshape(vp(:,ii),sz,sz)')
    axis off
    colormap hot

    hold on

    sg = up(:,ii);
    sg = normalize(sg,'range')*sz;
    %sg = filtfilt(b,a,sg);
    sg = resample(sg,sz,1024);
    plot(sg,'Color','green');
    
    hold on
    text(5,15,sprintf("#%d",ii));
    
    
    hold off
end

% filter 
figure;
nr = 6;nc=6;
tiledlayout(nr,nc,'TileSpacing','tight','Padding','none');

for ii=1:nr*nc
    nexttile
    imagesc(reshape(vpf(:,ii),sz,sz)')
    axis off
    colormap hot

    hold on

    sg = upf(:,ii);
    sg = normalize(sg,'range')*sz;
    %sg = filtfilt(b,a,sg);
    sg = resample(sg,sz,1024);
    plot(sg,'Color','green');
    
    hold on
    text(5,15,sprintf("#%d",ii));
    
    
    hold off
end

%% proof 2: some leaking harmonics
figure;
tt = tiledlayout(1,3);

faxis = (-512:511)*Fs/1024;
svaxis = 1:1024;

nexttile
ff = fftshift(1:1024);
ftdat = mag2db(abs(fft(up)));
imagesc(svaxis,faxis/1e9,ftdat(ff,:));hold on
xline([25 200 300],'LineStyle','--','HandleVisibility','off','LineWidth',2,'Alpha',0.2,'Color','black');hold off
colormap jet
hcb=colorbar;
hcb.Title.String = "(dB)";
set(gca,'FontSize',16,'FontName','Times New Roman')
%xlabel("# Number of singular values",'FontSize',16,'FontName','Times New Roman')
ylabel("Frequency (GHz)",'FontSize',16,'FontName','Times New Roman')
axis square
title("$|\mathcal{F} (U_{P})|$",'Interpreter', 'LaTeX')

nexttile
ff = fftshift(1:1024);
ftdat = mag2db(abs(fft(uh)));
imagesc(svaxis,faxis/1e9,ftdat(ff,:));
colormap jet
axis square
hcb=colorbar;
hcb.Title.String = "(dB)";
set(gca,'FontSize',16,'FontName','Times New Roman')
%xlabel("# Number of singular values",'FontSize',16,'FontName','Times New Roman')
ylabel("Frequency (GHz)",'FontSize',16,'FontName','Times New Roman')
title("$|\mathcal{F} (U_{\widetilde{H}})|$",'Interpreter', 'LaTeX')

% proof 3: similarity in sigma
nexttile
semilogy(diag(sp),'DisplayName','$\Sigma_P$','LineWidth',2);hold on
semilogy(diag(sh),'DisplayName','$\Sigma_{\widetilde{H}}$','LineWidth',2);hold on
xline([25 200 300],'LineStyle','--','HandleVisibility','off','LineWidth',2,'Alpha',0.2);hold off

axis square tight
legend('Interpreter','latex')
set(gca,'FontName','Times','FontSize',16)
ylabel("Energy",'FontSize',16,'FontName','Times New Roman')

xlabel(tt,"# Number of singular values",'FontSize',16,'FontName','Times New Roman')

%% Filter U
%[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
% [b,a] = butter(4,[100e6] / (Fs/2),'low');
% 
% upf = filtfilt(b,a,up);
% 
% space2Df = upf*(sp)*vp';


% figure
% tiledlayout(1,6);
% for ii=1:6
%     nexttile
%     plot(up(:,ii)); hold on
%     plot(upbpf(:,ii)); hold off
%     axis square tight
% end

%%
% figure
% imagesc(squeeze(max(abs(hilbert(space3D)))));
% axis image;colormap hot
% 
% nsphere = 9;
% onesphererad = 10;
% area = nsphere*2*pi*onesphererad^2;
% area*100 / 100^2
%% cmode
figure
imagesc(squeeze(max(abs(hilbert(space3D)))));
axis image;colormap jet
%% SNR
strong_signal = squeeze(space3D(:,155,200));
weak_signal = squeeze(space3D(:,9,28));
ns = squeeze(space3D(:,44,110));
faxis = (1:512)*Fs/1024;
figure;
ft = abs(fft(strong_signal));
plot(faxis,mag2db(ft(1:512)));
hold on
ft = abs(fft(weak_signal));
plot(faxis,mag2db(ft(1:512)));
hold on
ft = abs(fft(ns));
plot(faxis,mag2db(ft(1:512)));
hold off

%% compare wU wV
range = @(x)(max(x) - min(x));

D = awgn(space3D,60);
% re-arrange
imsize = size(D,2);
dat2d = reshape(D,1024,imsize^2);
[u,s,v] = svd(dat2d,"econ");

[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
u = filtfilt(b,a,u);

% wU
uu = normalize(u,1,"range")-0.5;
uu = sum(uu.^2);
wU = normalize(-uu,'range');

% wV
v3d = reshape(v,imsize,imsize,1024);
ss = pagesvd(v3d);
vv = sum(squeeze(ss));
wV = normalize(-vv,'range');

%% compare
cmodefunc = @(x)(squeeze(max(abs(hilbert(x)))));
Pu = u*(s.*diag(wU))*v';
Pv = u*(s.*diag(wU))*v';
Pbpf = u*s*v';

Pu3d = reshape(Pu,1024,sz,sz);
Pv3d = reshape(Pv,1024,sz,sz);
Pbpf3d = reshape(Pbpf,1024,sz,sz);

figure;
tiledlayout(1,4)
nexttile
imagesc(cmodefunc(D));axis image
nexttile
imagesc(cmodefunc(Pbpf3d));axis image
nexttile
imagesc(cmodefunc(Pu3d));axis image
nexttile
imagesc(cmodefunc(Pv3d));axis image


%% SNR improvement
dat = Pbpf3d;
strong_signal = squeeze(dat(:,82,83));
weak_signal = squeeze(dat(:,80,11));
ns = squeeze(dat(:,48,21));
faxis = (1:512)*Fs/1024;

figure;
ft = abs(fft(strong_signal));
plot(faxis,mag2db(ft(1:512)));
hold on
ft = abs(fft(weak_signal));
plot(faxis,mag2db(ft(1:512)));
hold on
ft = abs(fft(ns));
plot(faxis,mag2db(ft(1:512)));
hold off

figure;
taxis = (1:1024)*1e6/Fs;
plot(taxis,strong_signal);
hold on
plot(taxis,weak_signal);
hold on
plot(taxis,ns);
hold off
axis tight