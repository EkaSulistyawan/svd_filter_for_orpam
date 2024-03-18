%% plot
titlelist = [['Raw'],"BPF","BPF+SVD wU", "BPF+SVD wV",...
    "Sym6-SURE","Sym4-DWT","Sym4-MODWT"];
figure("Position",[2012,113,1623,916])
tl = tiledlayout(4,7,'TileSpacing','tight');

seldat.raw = space3D;
seldat.bpf = space3Dbpf;
seldat.wu = denoisedSVDwU;
seldat.wv = denoisedSVDwV;
seldat.sym6 = denoisedPaper1Wavelet;
seldat.sym4dwt = denoisedPaper2DWT;
seldat.sym4modt = denoisedPaper2MODWT;

varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    hold on
    % view roi
    % xbox = [sgx(1), sgx(end), sgx(end), sgx(1)  , sgx(1)];
    % ybox = [sgy(1), sgy(1)  , sgy(end), sgy(end), sgy(1)];
    % plot(xbox,ybox,'Color','black','LineWidth',1.5); hold on
    % xbox = [bgx(1), bgx(end), bgx(end), bgx(1)  , bgx(1)];
    % ybox = [bgy(1), bgy(1)  , bgy(end), bgy(end), bgy(1)];
    % plot(xbox,ybox,'Color','green','LineWidth',1.5); hold off
    line([1 120],[71 71],'LineWidth',2,'Color','g','LineStyle','--')
    % line([1 33],[73 73],'LineWidth',2,'Color','c','LineStyle','--')
    line([56 76],[14 14],'LineWidth',2,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+14)
xaxis = (1:33)*250 / 1000;
plot(xaxis,im(1:33,71),'LineWidth',2)
axis tight
ylim([0.05 0.8])
xlabel(['Length (' char(181) 'm)'])
xlabel('Length ($\mu$m)','Interpreter','latex');
ylabel(' Norm. Intensity (a.u.)')
set(gca,'FontName','Times New Roman','FontSize',12)


nexttile(i+7)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:120)*250 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
% xlabel(['Length (' char(181) 'm)'])
% ylabel(['Depth (' char(181) 'm)'])
xlabel('Length ($\mu$m)','Interpreter','latex');
ylabel('Depth ($\mu$m)','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+21)
taxis = (1:1024)*1e6/Fs;
dat3d = seldat.(varnames{i});
plot(taxis,dat3d(:,9,55),'DisplayName','Signal','Color','Red','LineWidth',2);hold on;
plot(taxis,dat3d(:,37,27),'DisplayName','Unwanted','Color','Black','LineWidth',2);hold off
xlabel("Time ($\mu$s)",'Interpreter','latex')
ylabel("Intensity (a.u.)")
ylim([-3e-4 5e-4]);
set(gca,'FontName','Times New Roman','FontSize',12)


end


% cbh = colorbar;
% cbh.Layout.Tile = 'north';
% put legend
leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
leg.Layout.Tile = 'south';
ax = gcf;
%exportgraphics(ax,'./exported_images/fig10_rbc_v2.emf')

%% only two line
titlelist = [['Raw'],"BPF","BPF+SVD wU", "BPF+SVD wV",...
    "Sym6-SURE","Sym4-DWT","Sym4-MODWT"];
figure("Position",[2012,113,1623,916])
tl = tiledlayout(2,7,'TileSpacing','tight');

seldat.raw = space3D;
seldat.bpf = space3Dbpf;
seldat.wu = denoisedSVDwU;
seldat.wv = denoisedSVDwV;
seldat.sym6 = denoisedPaper1Wavelet;
seldat.sym4dwt = denoisedPaper2DWT;
seldat.sym4modt = denoisedPaper2MODWT;

varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    hold on
    % view roi
    xbox = [sgx(1), sgx(end), sgx(end), sgx(1)  , sgx(1)];
    ybox = [sgy(1), sgy(1)  , sgy(end), sgy(end), sgy(1)];
    plot(xbox,ybox,'Color','black','LineWidth',1.5); hold on
    xbox = [bgx(1), bgx(end), bgx(end), bgx(1)  , bgx(1)];
    ybox = [bgy(1), bgy(1)  , bgy(end), bgy(end), bgy(1)];
    plot(xbox,ybox,'Color','green','LineWidth',1.5); hold off
    line([1 200],[9 9],'LineWidth',2,'Color','g','LineStyle','--')
    % line([68 165],[9 9],'LineWidth',2,'Color','c','LineStyle','--')
    line([10 30],[190 190],'LineWidth',2,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)


nexttile(i+7)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:120)*250 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)

end

%% ROI
figure;
i = 1;
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');colormap hot;axis off
if(i == 1)
    hold on
    %view roi
    xbox = [sgx(1), sgx(end), sgx(end), sgx(1)  , sgx(1)];
    ybox = [sgy(1), sgy(1)  , sgy(end), sgy(end), sgy(1)];
    plot(xbox,ybox,'Color','black','LineWidth',1.5); hold on
    xbox = [bgx(1), bgx(end), bgx(end), bgx(1)  , bgx(1)];
    ybox = [bgy(1), bgy(1)  , bgy(end), bgy(end), bgy(1)];
    plot(xbox,ybox,'Color','green','LineWidth',1.5); hold off
    %line([1 120],[71 71],'LineWidth',2,'Color','g','LineStyle','--')
    % line([1 33],[73 73],'LineWidth',2,'Color','c','LineStyle','--')
    %line([56 76],[14 14],'LineWidth',2,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)
axis image;
%% show image collection
load("rbc_denoising.mat")
load("rbc_denoising_USV.mat",'u','v','s','ubpf','vbpf','sbpf','sz')
%
up = u; vp = v; sp = s;

% recalculate weight
[wu,wv,ws] = getweight(up,vp,sp);

%% no filter 
figure('Position',[59,65,1061,914]);
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
    plot(sg,'Color','Cyan','LineWidth',2);
    
    hold on
    bgc = [1 1 1 0.75];
    text(5,15,sprintf("#%d",ii),'FontSize',18,'BackgroundColor',bgc);
    text(65,15,sprintf("$\\Sigma=%.2f$",ws(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,85,sprintf("$w_U=%.2f$",wu(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,110,sprintf("$w_V=%.2f$",wv(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    
    
    hold off
end

% filter 
upf = ubpf; vpf = vbpf; spf = sbpf;
[wu,wv,ws] = getweight(upf,vpf,spf);

figure('Position',[59,65,1061,914]);
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
    plot(sg,'Color','Cyan','LineWidth',2);
    
    hold on
    bgc = [1 1 1 0.75];
    text(5,15,sprintf("#%d",ii),'FontSize',18,'BackgroundColor',bgc);
    text(65,15,sprintf("$\\Sigma=%.2f$",ws(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,85,sprintf("$w_U=%.2f$",wu(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,110,sprintf("$w_V=%.2f$",wv(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    
    
    hold off
end

%%
i=4;
figure;
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');colormap hot;axis square off

figure
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:120)*250 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
axis off
% xlabel(['Length (' char(181) 'm)'])
% ylabel(['Depth (' char(181) 'm)'])
% xlabel('Length ($\mu$m)','Interpreter','latex');
% ylabel('Depth ($\mu$m)','Interpreter','latex');
% set(gca,'FontName','Times New Roman','FontSize',12)
%% for graphical abstract
im = cmodefunc(space3Dbpf);
sg = space3Dbpf(:,10,74);
sg = normalize(sg,'range')*sz;
%sg = filtfilt(b,a,sg);
sg = resample(sg,sz,1024);

figure;
imagesc(im');hold on
plot(sg,'Color','Cyan','LineWidth',4);hold off
colormap hot
axis image off

%%
im = cmodefunc(denoisedSVDwV);
sg = denoisedSVDwV(:,10,74);
sg = normalize(sg,'range')*sz;
%sg = filtfilt(b,a,sg);
sg = resample(sg,sz,1024);

figure;
imagesc(im');hold on
plot(sg,'Color','Cyan','LineWidth',4);hold off
colormap hot
axis image off

%%
im = cmodefunc(denoisedSVDwU);
sg = denoisedSVDwU(:,10,74);
sg = normalize(sg,'range')*sz;
%sg = filtfilt(b,a,sg);
sg = resample(sg,sz,1024);

figure;
imagesc(im');hold on
plot(sg,'Color','Cyan','LineWidth',4);hold off
colormap hot
axis image off

%% graphical abstract v2
% raw
im = cmodefunc(space3D);
sg = space3D(:,10,74);
sg = normalize(sg,'range')*sz;
%sg = filtfilt(b,a,sg);
sg = resample(sg,sz,1024);
figure;
imagesc(im');hold on
plot(sg,'Color','Cyan','LineWidth',4);hold off
colormap hot
axis image off


isel = [1,4,6];
for i = 1:3
    figure
    ii = isel(i);
    imagesc(reshape(vpf(:,ii),sz,sz)')
    axis off square
    colormap hot
    
    hold on
    
    sg = upf(:,ii);
    sg = normalize(sg,'range')*sz;
    %sg = filtfilt(b,a,sg);
    sg = resample(sg,sz,1024);
    plot(sg,'Color','Cyan','LineWidth',6);
end

isel = [1,4,6];
for i = 1:3
    figure
    ii = isel(i);
    imagesc(reshape(vpf(:,ii),sz,sz)')
    axis off square
    colormap hot
    
    hold on
    
    sg = upf(:,ii);
    sg = normalize(sg,'range')*sz;
    mm = mean(sg);
    sg2 = normalize(sg,'range')*sz*wu(ii);
    mm2 = mean(sg2);
    %sg = filtfilt(b,a,sg);
    sg = resample(sg2+(mm-mm2),sz,1024);
    sg(1) = mean(sg(2:5));
    plot(sg,'Color','Cyan','LineWidth',6);
end
%% SNR calculation
% a = denoisedSVDwV(:,sgx(1),sgy(1));
% b = denoisedSVDwV(:,bgx(1),bgy(1));
a = space3Dbpf(:,sgx(1),sgy(1));
b = space3Dbpf(:,bgx(1),bgy(1));
lowlim = 10e6;
hghlim = 100e6;
faxis = (1:512)*5e9/1024;
idx = find((faxis > lowlim) & (faxis < hghlim));
fta = abs(fft(a));
ftb = abs(fft(b));
maxdbsg = max(fta(idx));
maxdbns = max(ftb(idx));

figure('Position',[680,113,407,765])
tiledlayout(2,1);
nexttile
taxis = (1:1024)*1e6/Fs;
plot(taxis,a,'Color','Black','LineWidth',lw,'DisplayName','Signal');hold on
plot(taxis,b,'Color','Red','LineWidth',lw,'DisplayName','Signal');hold off
axis  tight
xlabel("Time (us)")
ylabel("Amplitude (a.u.)")
set(gca,'FontName','Times New Roman','FontSize',20)

nexttile
lw = 2;
plot(faxis*1e-9,mag2db(fta(1:512)),'Color','Black','LineWidth',lw,'DisplayName','Signal');hold on
plot(faxis*1e-9,mag2db(ftb(1:512)),'Color','Red','LineWidth',lw,'DisplayName','Background');hold on
yline(mag2db(maxdbsg),'Color',[0.1 0.1 0.1 0.9],'LineWidth',lw,'DisplayName','Signal level'); hold on
yline(mag2db(maxdbns),'Color',[0.5 0.1 0.1 0.9],'LineWidth',lw,'DisplayName','Background level'); hold off

text(faxis(200)*1e-9,(mag2db(maxdbsg)+mag2db(maxdbns))/2,sprintf("SNR %.2f dB",mag2db(maxdbsg)-mag2db(maxdbns)),...
    'FontName','Times New Roman','FontSize',15)

xlabel("Frequency (GHz)")
ylabel("Magnitude (dB)")
set(gca,'FontName','Times New Roman','FontSize',20)

axis tight
ylim([-100 -20])
%%
function p = cmodefunc(x)
    p = squeeze(max(abs(hilbert(x))));
    p = (p - min(p(:))) / range(p(:));
    p(p > 0.45) = 0.45;
    p = (p - min(p(:))) / range(p(:));
end

function p = bmodefunc(x)
    hb = abs(hilbert(x));
    p = squeeze(hb(:,:,71));
    p = (p - min(p(:))) / range(p(:));
end

function [wU,wV,wS] = getweight(u,v,s)
    imsize = round(sqrt(size(v,1)));

    uu = normalize(u,1,"range")-0.5;
    uu = sum(uu.^2);
    wU = normalize(-uu,'range');
    %wU = wU' .* diag(s);
    
    % wV
    v3d = reshape(v,imsize,imsize,1024);
    ss = pagesvd(v3d,'econ');
    vv = sum(squeeze(ss));
    wV = normalize(-vv,'range');
    %wV = wV' .* diag(s);

    wS = diag(s);
end

