%% synthetic noise
dat = load("C:\Users\Eka\Documents\SharingPoint\data\cellular\compiled\04092023_USAF.mat").space3D;
dat = dat(:,1:size(dat,2)-1,1:size(dat,2)-1);
imsz = size(dat,2);

cmode = squeeze(max(abs(hilbert(dat))));
figure;imagesc((cmode));axis image off;colormap hot

% 45 ~ 1 dB SNR
% 60 ~ 10 dB SNR
datns = awgn(dat,60);
valsnr = get_snr(datns,'Flat');
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
space3Dbpf = filtfilt(b,a,datns);

dat2d = reshape(space3Dbpf,1024,imsz*imsz);


%% Other analysis
dat = load("C:\Users\Eka\Documents\SharingPoint\data\cellular\compiled\23022023_rbc.mat").space3D;
dat = dat(:,1:size(dat,2)-1,1:size(dat,2)-1);

cmode = squeeze(max(abs(hilbert(dat))));
figure;imagesc((cmode));axis image off;colormap hot
%figure;imagesc((cmode));axis image off;colormap hot
%%
bmode = (abs(hilbert(dat)));
bmode = bmode(:,138,:);
figure;imagesc(squeeze(bmode));axis off;colormap hot
%%
figure('Position',[680,458,838,420])
view3Dpa(abs(hilbert(dat)),0.4,10);
view([44.182274835291560,58.087548541656105])
colormap hot
xlim([1 Inf])

%% tear down
dat2d = reshape(datns,1024,100*100);
[u,s,v] = svd(dat2d,'econ');


%% on terrible measurement
dat = load("C:\Users\Eka\Documents\SharingPoint\data\cellular\compiled\21032023_rbc_noavg.mat").space3D;
dat = dat(:,1:size(dat,2)-1,1:size(dat,2)-1);
ref = load("C:\Users\Eka\Documents\SharingPoint\data\cellular\compiled\21032023_rbc_ref.mat").space3D;

ref = ref(:,1:size(ref,2)-1,1:size(ref,2)-1);
for ii=39:41
    ref(:,ii,46) = ref(:,ii,50);
end


Fs = 5e9;
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
space3Dbpf = filtfilt(b,a,dat);

Fs = 5e9;
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
refbpf = filtfilt(b,a,ref);

%%
cmode = squeeze(max(abs(hilbert(refbpf))));
figure;imagesc((cmode));axis image off;colormap hot

%% synthetic noise
close all

refns = awgn(ref,65);
Fs = 5e9;
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
refnsbpf = filtfilt(b,a,refns);
%
[snrr0ori,cnrr0ori] = get_metrics(refbpf,'RBC');
%%
[snrr1bpf,cnrr1bpf] = get_metrics(refnsbpf,'RBC');
cmodebpf = squeeze(max(abs(hilbert(refnsbpf))));
figure;imagesc((cmodebpf));axis image off;colormap hot

figure(Position=[651,689,560,214])
plot(refns(:,137,116),DisplayName='wo BPF',Color=[.5, .5, .5],LineWidth=2);hold on
plot(refnsbpf(:,137,116),DisplayName='w BPF',Color='r',LineWidth=2);hold off
axis tight
% xlabel('Time (ns)')
% ylabel('Norm. Intensity (a.u.)')
ylim([-2e-3, 2e-3])
xticklabels([])
yticklabels([])

%
paper1emd = paper_denoising(refnsbpf, "paper-1-emd-mi");
[snrr2emd,cnrr2emd] = get_metrics(paper1emd,'RBC');
cmodeemd = squeeze(max(abs(hilbert(paper1emd))));
figure;imagesc((cmodeemd));axis image off;colormap hot

%
paper2mowv = paper_denoising(refnsbpf, "paper-2-modwt");
[snrr3modwt,cnrr3modwt] = get_metrics(paper2mowv,'RBC');
cmodemodwt = squeeze(max(abs(hilbert(paper2mowv))));
figure;imagesc((cmodemodwt));axis image off;colormap hot

%
[ours,~,~,weight] = svd_denoising_ius(refnsbpf,0);
[snrr4ours,cnrr4ours] = get_metrics(ours,'RBC');
cmodeours = squeeze(max(abs(hilbert(ours))));
figure;imagesc((cmodeours));axis image off;colormap hot

%%
[snrr5bpfnoavg,cnrr5bpfnoavg] = get_metrics(space3Dbpf,'RBC');

%% test paper denoising
% sig = dat(:,137,116);
% sigbpf = space3Dbpf(:,137,116);
paper1emd = paper_denoising(dat, "paper-1-emd-mi");
% paper1wv = paper_denoising(dat, "paper-1-wavelet");
% paper2wv = paper_denoising(dat, "paper-2-dwt");
paper2mowv = paper_denoising(dat, "paper-2-modwt");
[ours,~,wvnorm,weight]=svd_denoising_ius(space3Dbpf);

%%
s = svd(reshape(space3Dbpf,1024,200*200),'econ','vector');
loglog(s,LineWidth=5);hold on
loglog(wvnorm,LineWidth=3,Color=[0, 1, 0, 0.8]);hold on
xline(3,linestyle='--',LineWidth=3);hold off
grid on
axis square
xticks([3])
ylim([1e-3, 1])
set(gca,'FontName','Aptos','FontSize',20)


%%
figure;imagesc(squeeze(max(abs(hilbert(space3Dbpf)))));axis image off;colormap hot

%%
[a,b] = get_metrics(space3Dbpf,'RBC');
%% signal at brightest RBC
close all
% figure(Position=[651,689,560,214])
figure(Position=[438,602,765,214])
plot(refbpf(:,137,116),DisplayName='Signal',Color=[0.0 0.5 0.0],LineWidth=2);hold on
plot(refbpf(:,65,31),DisplayName='Noise',Color=[1.0 0.0 0.5],LineWidth=2);hold on
axis tight
% xlabel('Time (ns)')
% ylabel('Norm. Intensity (a.u.)')
ylim([-1e-3, 1e-3])
xticklabels([])
yticklabels([])
set(gca,'FontName','Aptos','FontSize',15)
% legend

% figure(Position=[651,689,560,214])
figure(Position=[438,602,765,214])
plot(space3Dbpf(:,137,116),DisplayName='Signal',Color=[0.0 0.5 0.0],LineWidth=2);hold on
plot(space3Dbpf(:,65,31),DisplayName='Noise',Color=[1.0 0.0 0.5],LineWidth=2);hold on
axis tight
% xlabel('Time (ns)')
% ylabel('Norm. Intensity (a.u.)')
ylim([-1e-3, 1e-3])
xticklabels([])
yticklabels([])
set(gca,'FontName','Aptos','FontSize',15)
% legend

% figure(Position=[651,689,560,214])
figure(Position=[438,602,765,214])
plot(paper1emd(:,137,116),DisplayName='Signal',Color=[0.0 0.5 0.0],LineWidth=2);hold on
plot(paper1emd(:,65,31),DisplayName='Noise',Color=[1.0 0.0 0.5],LineWidth=2);hold on
axis tight
% xlabel('Time (ns)')
% ylabel('Norm. Intensity (a.u.)')
ylim([-1e-3, 1e-3])
xticklabels([])
yticklabels([])
set(gca,'FontName','Aptos','FontSize',15)
% legend

% figure(Position=[651,689,560,214])
figure(Position=[438,602,765,214])
plot(paper2mowv(:,137,116),DisplayName='Signal',Color=[0.0 0.5 0.0],LineWidth=2);hold on
plot(paper2mowv(:,65,31),DisplayName='Noise',Color=[1.0 0.0 0.5],LineWidth=2);hold on
axis tight
% xlabel('Time (ns)')
% ylabel('Norm. Intensity (a.u.)')
ylim([-1e-3, 1e-3])
xticklabels([])
yticklabels([])
set(gca,'FontName','Aptos','FontSize',15)
% legend

% figure(Position=[651,689,560,214])
figure(Position=[438,602,765,214])
plot(ours(:,137,116),DisplayName='Signal',Color=[0.0 0.5 0.0],LineWidth=2);hold on
plot(ours(:,65,31),DisplayName='Noise',Color=[1.0 0.0 0.5],LineWidth=2);hold on
axis tight
% xlabel('Time (ns)')
% ylabel('Norm. Intensity (a.u.)')
ylim([-1e-3, 1e-3])
xticklabels([])
yticklabels([])
set(gca,'FontName','Aptos','FontSize',15)

%%

figure
imagesc(squeeze(max(abs(hilbert(paper2wv)))))
axis square off
colormap hot
%%
ttime = (1:1024) * 1e9 / Fs;
imagesc(reshape(ours,1024,200*200))
colormap gray
xticklabels([])
yticklabels([])
%% by paper

% plot(A0,'DisplayName','Sym6-SURE');hold on
ttime = (1:1024) * 1e9 / Fs;

figure(Position=[651,689,560,214])
plot(ttime,paper1emd,'DisplayName','EMD-MI',LineWidth=2);hold on
%plot(ttime,paper1wv,'DisplayName','[1]',LineWidth=2);hold on
% plot(ttime,paper2wv,'DisplayName','[2]',LineWidth=2);hold on
plot(ttime,paper2mowv,'DisplayName','MODWT',LineWidth=2);hold off
% plot(ttime,sigbpf,LineStyle='--',Color='r',DisplayName='Reference',LineWidth=1);hold off
axis tight
ylim([-2e-3, 2e-3])
% xlabel('Time (ns)')
% ylabel('Norm. Intensity (a.u.)')
xticklabels([])
yticklabels([])
set(gca,'FontName','Aptos','FontSize',20)
% legend(Orientation="horizontal",Location='southoutside')
%% by svd

[u,s,v] = svd(reshape(space3Dbpf,1024,200*200),'econ');
v3d = reshape(v',1024,200,200);

%%
loglog(s,LineWidth=4)
xline(4,LineStyle='--',LineWidth=1.2)
xline(40,LineStyle='--',LineWidth=1.2)
xticks([4, 40])
%yticks([1e-20, 1e-10, 1])
ylim([1e-3, 1])
set(gca,'FontName','Aptos','FontSize',20)
axis square
grid on

%%
close all
idx = 3;

vim = squeeze(v3d(idx,:,:));

figure;
imagesc(vim')
hold on

% get the signal
sg = u(:,idx);
sg = resample(sg,200,1024);
sg = normalize(sg,'range') * 200;
plot(sg,Color='cyan',LineWidth=5)
colormap hot
axis image off

sim = svd(vim,"econ","vector");
sim = sim';
%%
figure
% Shade the area below the curve
x = 1:length(sim);          % X-axis values
y = sim;                    % Y-axis values
fill([x fliplr(x)], [y zeros(1,length(y))], ones(1,3)*0.5, 'FaceAlpha', 0.8, 'EdgeColor', 'none')
hold on
plot(sim,Color='r',LineWidth=5)
hold off



% Set axis properties
xlim([1,200])
xticklabels([])
yticklabels([])
ylim([0, 0.15])

axis square
%%
figure('Position',[680,458,838,420])
view3Dpa(space3Dbpf,0.4,10);
view([44.182274835291560,58.087548541656105])
colormap hot
zlim([100 150])
xlim([1 Inf])

%%
[~,~,wv,weight]=svd_denoising_ius(space3Dbpf);
semilogx(weight)
%%
cmode = squeeze(max(abs(hilbert(refbpf))));
figure;imagesc((cmode));axis image off;colormap hot
%figure;imagesc((cmode));axis image off;colormap hot

figure('Position',[680,458,838,420])
view3Dpa(refbpf,0.4,10);
view([44.182274835291560,58.087548541656105])
colormap hot
zlim([100 150])
xlim([1 Inf])

%%
sz = 200;
Fs = 5e9;
kz = (1:1024)*1e9/Fs;
kx = (1:200)*20/200;

sgperi = mean(squeeze(refbpf(:,109:113,1)),2);
sgcent = mean(squeeze(refbpf(:,123:127,1)),2);

bmode = (abs(hilbert(refbpf)));
bmode = bmode(:,:,1);
figure;
imagesc(kx,kz*1e-9*1e6*1500,squeeze(bmode));
colormap hot
%ylim([100 150])
set(gca,'FontName','Times','FontWeight','bold','FontSize',20)
colorbar

figure;
imagesc(bmode);

figure;
plot(sgperi);hold on
plot(abs(hilbert(sgperi)),'Color',[1 0 1]);hold off
axis tight off


figure;
plot(sgcent);hold on
plot(abs(hilbert(sgcent)),'Color',[1 0 1]);hold off
axis tight off

%% compare
figure
plot(kz*1e-9*1e6*1500,abs(hilbert(sgcent)),'DisplayName','Center','LineWidth',2);hold on
plot(kz*1e-9*1e6*1500,abs(hilbert(sgperi)),'DisplayName','Periphery','LineWidth',2);hold on
xline(165,'LineWidth',2,'LineStyle','--','HandleVisibility','off');hold on
xline(193.8,'LineWidth',2,'LineStyle','--','HandleVisibility','off');hold on
xline(198.3,'LineWidth',2,'LineStyle','--','HandleVisibility','off');hold on
yline(2e-4,'LineWidth',2,'LineStyle','--','HandleVisibility','off')
axis tight
legend('Location','northwest')
set(gca,'FontName','Times','FontWeight','bold','FontSize',20)

plot(kz*1e-9*1e6*1500,((sgcent)),'DisplayName','Center','LineWidth',2);hold on
plot(kz*1e-9*1e6*1500,((sgperi)),'DisplayName','Periphery','LineWidth',2);hold off
axis tight
legend('Location','northwest')
set(gca,'FontName','Times','FontWeight','bold','FontSize',20)


%% fwhm
gg = fit(kz'*1e-9*1e6*1500,abs(hilbert(sgcent)),'gauss1');
sigma = gg.c1/sqrt(2); 
FWHM = 2*sqrt(2*log(2))*sigma;
fprintf("FWHM: %.2f um \n ",FWHM)

figure;
plot(gg); hold on

gg = fit(kz'*1e-9*1e6*1500,abs(hilbert(sgperi)),'gauss1');
sigma = gg.c1/sqrt(2); 
FWHM = 2*sqrt(2*log(2))*sigma;
fprintf("FWHM: %.2f um \n ",FWHM)
plot(gg);hold off
axis tight\
%% apply filter
tic
[res,wU,wV] = svd_denoising_ius(space3Dbpf,0);
toc
%% get SNR
[snrr,cnrr] = get_metrics(res,'RBC');
%
[snrrraw,cnrrraw] = get_metrics(space3Dbpf,'RBC');
%
[snrrref,cnrrref] = get_metrics(refbpf,'RBC');

%% 
cmoderef = squeeze(max(abs(hilbert(refbpf))));
cmoderes = squeeze(max(abs(hilbert(res))));
cmoderaw = squeeze(max(abs(hilbert(space3Dbpf))));
%%
figure;
imagesc(rescale(cmoderaw));colorbar

figure;
imagesc(rescale(cmoderes));colorbar

figure;
imagesc(rescale(cmoderef));colorbar


psnr(rescale(cmoderes),rescale(cmoderef))
psnr(rescale(cmoderaw),rescale(cmoderef))
%%
dat2d = reshape(ref,1024,200*200);
[u,s,v] = svd(dat2d,'econ');
v2d = reshape(v,200,200,1024);
 
%%
close all
for ii=1:5
%ii = 1;

figure;
imagesc(-(squeeze(v2d(:,:,ii))));axis image off;colormap hot

[u1,s1,v1] = svd((squeeze(v2d(:,:,ii))),'econ');
figure;
plot(diag(s1),'Color','Red','LineWidth',5); hold on
area([1:200],diag(s1),'FaceColor',0.5*[1 1 1],'EdgeColor','none');hold off
set(gca,'XTickLabels',[],'YTickLabels',[])
axis square off


figure;
plot(u(:,ii),'Color','Cyan','LineWidth',5);axis off
end
%%
close all
for ii=1:10
%ii = 1;

figure;
subplot(121);imagesc((squeeze(v2d(:,:,ii))));axis image off;colormap jet;
subplot(122);plot(filtfilt(b,a,u(:,ii)),'Color','Cyan','LineWidth',5);axis off
end
%% center periphery by svd
ii = 1;
sgcent = (squeeze(v2d(127,3,ii)))*u(:,ii);
sgperi = (squeeze(v2d(139,18,ii)))*u(:,ii);

figure
plot(kz*1e-9*1e6*1500,abs(hilbert(sgcent)),'DisplayName','Center','LineWidth',2);hold on
plot(kz*1e-9*1e6*1500,abs(hilbert(sgperi)),'DisplayName','Periphery','LineWidth',2);hold on
axis tight
legend('Location','northwest')
set(gca,'FontName','Times','FontWeight','bold','FontSize',20)

plot(kz*1e-9*1e6*1500,((sgcent)),'DisplayName','Center','LineWidth',2);hold on
plot(kz*1e-9*1e6*1500,((sgperi)),'DisplayName','Periphery','LineWidth',2);hold off
axis tight
legend('Location','northwest')
set(gca,'FontName','Times','FontWeight','bold','FontSize',20)
%%
cmode = squeeze(max(abs(hilbert(res))));
figure;imagesc((cmode));axis image off;colormap hot

figure('Position',[680,458,838,420])
view3Dpa(res,0.4,10);
view([44.182274835291560,58.087548541656105])
colormap hot
zlim([150 250])
xlim([1 Inf])

%% 
function status = view3Dpa(seldat,ths,ds)
    sz = size(seldat,2);
    dat2d = reshape(seldat,1024,sz*sz);
    dat2d = abs(hilbert(dat2d));
    dat1d = dat2d(:);
    
    dat1d = normalize(dat1d,'range');
    
    k = find(dat1d > ths);
    
    Fs = 5e9;
    kz = mod(k,1024)*1e9/Fs;
    kx = mod(ceil(k / 1024),sz)*200/1000;
    ky = ceil(ceil(k/1024)/sz)*200/1000;
    
    S = dat1d(k);
    C = dat1d(k);
    
    scatter3(kx(1:ds:end),ky(1:ds:end),kz(1:ds:end),S(1:ds:end).*10,C(1:ds:end));colormap hot
    colorbar

    status = true;
end

function [a,b] = get_metrics(D,datused)
    global_param;
    if datused == "Flat"
        sgx = sgx_flat2; 
        sgy = sgy_flat2;
        bgx = bgx_flat2; 
        bgy = bgy_flat2;
    elseif datused == "RBC"
        sgx = sgx_rbcIUS; 
        sgy = sgy_rbcIUS;
        bgx = bgx_rbcIUS; 
        bgy = bgy_rbcIUS;
    end
    
    sg = reshape(D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
    bg = reshape(D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
    
    taxis = (1:1024)*1e9/5e9;

    figure(Position=[438,602,765,214])
    alph = 0.3;
    plot_shade(taxis,squeeze(sg),'Signal',[0.0 0.5 0.0],alph);hold on
    plot_shade(taxis,squeeze(bg),'Noise',[1.0 0.0 0.5],alph);hold off
    axis tight
    set(gca,'FontName','Times','FontSize',25,'FontWeight','Bold')
    % ylim([-6 6]) % e-4
    ylim([-1 1]*10)
    %xticklabels([])
    yticklabels([])


    a = snr2d(sg,bg);
    b = get_cnr(D,sgx,sgy,bgx,bgy);
end

function [outp] = plot_shade(x,dat,nm,c1,ca)
    avg = mean(squeeze(dat),2)*1e4;
    stdv = std(squeeze(dat),[],2)*1e4;

    curve1 = avg + stdv;
    curve2 = avg - stdv;

    x2 = [x, fliplr(x)];
    inBetween = [curve1; flipud(curve2)];
    plot(x, avg,'Color', c1, 'LineWidth', 4,'DisplayName',strcat(['Avg.' nm]));hold on
    fill(x2, inBetween,c1,'FaceAlpha',ca,'EdgeColor','none','DisplayName',strcat(['Std.' nm]));hold on
    
end

function [outp] = snr2d(a,b)
    sz = size(a,2);
    outp = 0;
    for i =1:sz
        outp = outp + custom_snr(a(:,i),b(:,i),i);
    end
    outp = outp / sz;
end

function [a] = range(p)
    a = max(p,[],"all") - min(p,[],'all');
end


function [outp] = custom_snr(a,b,i)
    lowlim = 10e6;
    hghlim = 100e6;

    faxis = (1:512)*5e9/1024;

    idx = find((faxis > lowlim) & (faxis < hghlim));
    fta = abs(fft(a));
    ftb = abs(fft(b)); 
    maxdbsg = max(fta(idx));
    maxdbns = std(ftb(idx));
    
    % if i==1
    %     disp(i)
    %     figure;
    %     plot(fta(idx));hold on
    %     plot(ftb(idx));hold off
    %     title(i)
    % end
    outp = mag2db(maxdbsg/maxdbns);
    % fprintf("%.2f\n",outp);
end

function [outp] = snrv2(a,b,i)
    lowlim = 10e6;
    hghlim = 100e6;

    faxis = (1:512)*5e9/1024;

    idx = find((faxis > lowlim) & (faxis < hghlim));
    fta = abs(fft(a));
    ftb = abs(fft(b)); 
    maxdbsg = max(fta(idx));
    maxdbns = std(ftb(idx));
    
    outp = mag2db(maxdbsg/maxdbns);

end

function [a] = get_cnr(D,sgx,sgy,bgx,bgy)

    hb = abs(hilbert(D));
    cmode = squeeze(max(hb));
    sgmu = mean(cmode(sgx,sgy),"all");
    bgmu = mean(cmode(bgx,bgy),"all");
    sgsd = std(cmode(sgx,sgy),0,"all");
    bgsd = std(cmode(bgx,bgy),0,"all");
    a = (sgmu - bgmu) / (sqrt((sgsd^2 + bgsd^2) / 2));
    a = mag2db(a);
end