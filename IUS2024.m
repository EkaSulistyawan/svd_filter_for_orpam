%% synthetic noise
dat = load("C:\Users\Eka\Documents\SharingPoint\data\cellular\compiled\04092023_USAF.mat").space3D;
dat = dat(:,1:size(dat,2)-1,1:size(dat,2)-1);

% 45 ~ 1 dB SNR
% 60 ~ 10 dB SNR
datns = awgn(dat,60);
valsnr = get_snr(datns,'Flat');
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
space3Dbpf = filtfilt(b,a,datns);

[res] = svd_denoising_ius(space3Dbpf,0);


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

%% signal at brightest RBC
sig = space3Dbpf(:,137,115);
plot(sig)
axis tight

%%
emd(sig)

%%
figure('Position',[680,458,838,420])
view3Dpa(space3Dbpf,0.4,10);
view([44.182274835291560,58.087548541656105])
colormap hot
zlim([100 150])
xlim([1 Inf])

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
axis tight
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

    figure
    alph = 0.3;
    plot_shade(taxis,squeeze(sg),'Signal',[0.0 0.5 0.0],alph);hold on
    plot_shade(taxis,squeeze(bg),'Noise',[1.0 0.0 0.5],alph);hold off
    axis tight
    set(gca,'FontName','Times','FontSize',25,'FontWeight','Bold')
    ylim([-6 6])


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