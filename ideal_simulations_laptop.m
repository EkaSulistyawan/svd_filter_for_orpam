clear all 
close all


% Simulation specs
N = 5000;
T = 1024;
l = 128;   % PA length
xaxis = 1:N;
Fs = 5e9;
taxis = (1:T)*1e6/Fs;

tiledlayout(1,4,'Padding','none');
% Construct h
nshape = linspace(0.5,-0.5,l);
h = [zeros(1,(T-l)/2),nshape,zeros(1,(T-l)/2)];

[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
% h = filtfilt(b,a,h);

nexttile
plot(taxis,h,'LineWidth',2,'Color','black');
axis square
xlabel("Time ($\mu$s)",'Interpreter','latex')
ylabel("Amplitude (a.u.)")
set(gca,'FontSize',16,'FontName','Times New Roman')
title("{\it h(\it{t})}")

% Construct H
H = zeros(T,T);
for i=1:T
    H(i,:) = circshift(h,i);
end
H = circshift(H,512,2);

nexttile
imagesc(taxis,taxis,H);
xlabel("Time ($\mu$s)",'Interpreter','latex')
ylabel("Time ($\mu$s)",'Interpreter','latex')
axis image
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it H")


% construct R
R = zeros(T,N);

% rule 1
kpctg = 1/T;
rnd = randperm(T*N,kpctg*T*N);
R(rnd) = 1;

% rule 2
rule2 = 1:N;
R(512,rule2) = 1;

% rule 3
% rule3 = randperm(N,1);
% rule3 = rule3:(rule3+0.04*N);
% R(200,rule3) = 1;

% rule 4
% rule4 = randperm(N,1);
% rule4 = rule4 :(rule4+0.08*N);
% R(800,rule4) = 1;

% rule 5
rule5 = 1000:2000;
R(400,rule5) = 1;

% rule 6
rule5 = 4000:4500;
R(800,rule5) = 1;

nexttile
imagesc(xaxis,taxis,R);
xlabel(sprintf("Lateral \nPosition (a.u.)"))
ylabel("Time ($\mu$s)",'Interpreter','latex')
xticklabels([1000 3000 5000])
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it R")

% construct P
P = H*R;

nexttile
imagesc(xaxis,taxis,P);
xlabel(sprintf("Lateral \nPosition (a.u.)"))
ylabel("Time ($\mu$s)",'Interpreter','latex')
xticklabels([1000 3000 5000])
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it P")
%% P only
figure
xaxis = 1:size(P,2);
taxis = (1:size(P,1))*1e6/fs;
imagesc(xaxis,taxis,P);
axis square
xlabel("Lateral Position (a.u.)")
ylabel("Time ($\mu$s)",'Interpreter','latex')
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it P")
%% compute SVD
[UP,SP,VP] = svd(P);
[UH,SH,VH] = svd(H);
[UR,SR,VR] = svd(R);

%% proof 1
fewranks = 4;
gt = UP(:,1:fewranks);
eq5a= H*UR(:,1:fewranks);


figure;
tiledlayout(1,fewranks)
for i=1:fewranks
    nexttile
    plot(eq5a(:,i),'LineStyle','-','LineWidth',1.5,'Color','Black','DisplayName','H U_R');hold on
    plot(gt(:,i),'LineStyle','-','LineWidth',1.5,'Color','Red','DisplayName','U_P');hold off
    axis square
    set(gca,'FontSize',16,'FontName','Times New Roman')
    title(sprintf("Singular Vector #%d",i),'FontSize',14)
end
%leg = legend('Orientation', 'Horizontal','FontSize',14);

%% show HUP and VR
fewranks = 4;
figure;
tt = tiledlayout(2,1);

ttt = tiledlayout(tt,1,4);
ttt.Layout.Tile = 1;
for i=1:fewranks
    nexttile(ttt)
    plot(taxis,gt(:,i),'LineStyle','-','LineWidth',1.5,'Color','Red','DisplayName','$U_P$');hold on
    plot(taxis,eq5a(:,i),'LineStyle','-','LineWidth',1.5,'Color','Black','DisplayName','$H U_R$');hold off
    axis square
    set(gca,'FontSize',16,'FontName','Times New Roman')
    title(sprintf("#%d",i),'FontSize',14)
end
xlabel(ttt,'Time ($\mu$s)','Interpreter','latex','FontSize',16,'FontName','Times New Roman')
lgd = legend('Orientation',"vertical",'Interpreter','latex');
lgd.Layout.Tile = 'east';

eq5b = R'*VH;
mav = ones(1,10);

ttt = tiledlayout(tt,1,4);
ttt.Layout.Tile = 2;
for sel =1:fewranks
    nexttile(ttt)
    %plot(conv(abs(eq5b(:,sel)),mav),'Color',[0 0 0 0.2],'DisplayName','$R^T V_H$','LineWidth',1.5);hold on
    plot(conv(abs(VP(:,sel)),mav),'Color',[1 0 0 0.7],'DisplayName','$V_P$','LineWidth',1.5);hold on
    plot(conv(abs(VR(:,sel)),mav),'Color',[0 0 0 0.7],'DisplayName','$V_R$','LineWidth',1.5);hold off
    
    set(gca,'FontSize',16,'FontName','Times New Roman')
    axis square tight
    ylim([0 1])
end
xlabel(ttt,'Lateral Position (a.u.)','FontSize',16,'FontName','Times New Roman')
lgd = legend('Orientation',"vertical",'Interpreter','latex');
lgd.Layout.Tile = 'east';

ylabel(tt,'Amplitude (a.u.)','FontSize',16,'FontName','Times New Roman')

%% Proof 2 relation UP and UH and sigma
Fs = 5e9;
figure;
tt = tiledlayout(1,3);
nexttile
faxis = (-512:511)*Fs/1024;
ff = fftshift(1:1024);
ft = mag2db(abs(fft(UP)));
imagesc(1:1024,faxis/1e9,ft(ff,:));hold on
colormap jet
axis square
hcb=colorbar;
hcb.Title.String = "(dB)";
clim([-300 0])
set(gca,'FontSize',16,'FontName','Times New Roman')
title("$|\mathcal{F} (U_P)|$",'Interpreter', 'LaTeX')
ylabel("Frequency (GHz)",'FontSize',16,'FontName','Times New Roman')


nexttile
faxis = (-512:511)*Fs/1024;
ff = fftshift(1:1024);
ft = mag2db(abs(fft(UH)));
imagesc(1:1024,faxis/1e9,ft(ff,:));hold off
colormap jet
axis square
hcb=colorbar;
hcb.Title.String = "(dB)";
clim([-300 0])
set(gca,'FontSize',16,'FontName','Times New Roman')
title("$|\mathcal{F} (U_H)|$",'Interpreter', 'LaTeX')
ylabel("Frequency (GHz)",'FontSize',16,'FontName','Times New Roman')

% sigma
nexttile
semilogy(diag(SP),'DisplayName','$\Sigma_P$','LineWidth',1.5);hold on;
semilogy(diag(SH),'DisplayName','$\Sigma_H$','LineWidth',1.5);hold on;
semilogy(diag(SR),'DisplayName','$\Sigma_R$','LineWidth',1.5);hold off;
legend('Interpreter','latex','Location','southeast')
set(gca,'FontSize',16,'FontName','Times New Roman')
% xlabel('Singular values','FontSize',16,'FontName','Times New Roman')
ylabel('Energy','FontSize',16,'FontName','Times New Roman')

axis tight square
xlabel(tt,"# Number of singular values",'FontSize',16,'FontName','Times New Roman')
% axes('position',[.3 .4 .3 .3])
% box on
% 
% tmp = diag(SP);semilogy(tmp(1:20),'DisplayName','$\Sigma_P$','LineWidth',1.2);hold on;
% tmp = diag(SH);semilogy(tmp(1:20),'DisplayName','$\Sigma_H$','LineWidth',1.2);hold on;
% tmp = diag(SR);semilogy(tmp(1:20),'DisplayName','$\Sigma_R$','LineWidth',1.2);hold off;
% set(gca,'FontSize',14,'FontName','Times New Roman')


%% effect of noise
fth = mag2db(abs(fft(h)));
signal_level = max(fth);

noised = awgn(h,10);
n = noised - h;
ftnoised = mag2db(abs(fft(noised)));
ftn = mag2db(abs(fft(n)));

figure;
plot(noised);hold on
plot(h);hold off

figure;
plot(ftnoised);hold on
plot(ftn);hold off
xlim([0 512])
title("Frequency")

% %% add noise
% nslv = 10;
% PNP = awgn(P,nslv);
% PNH = awgn(H,nslv)*R;
% PNR = H*awgn(R,nslv);
% Ploc = P;
% 
% bpfon = true;
% if bpfon
%     Fs = 5e9;
%     [b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
%     PNR = filtfilt(b,a,PNR);
%     PNP = filtfilt(b,a,PNP);
%     PNH = filtfilt(b,a,PNH);
%     Ploc = filtfilt(b,a,P);
% end
% 
% [UNP,SNP,VNP] = svd(PNP);
% [UNR,SNR,VNR] = svd(PNR);
% [UNH,SNH,VNH] = svd(PNH);
% [UP,SP,VP] = svd(Ploc);
% 
% figure;
% tiledlayout(1,2)
% nexttile
% plot(UP(:,1:5));axis tight square
% nexttile
% plot(UNP(:,1:5));axis tight square
% % nexttile
% % plot(UNH(:,1:5));axis tight square
% % nexttile
% % plot(UNR(:,1:5));axis tight square
% 
% fs = 5e9;
% [b,a] = butter(4,[10e6 100e6] / (fs/2),'bandpass');
% 
% figure;
% tiledlayout(1,2)
% nexttile
% imagesc(P);axis square
% nexttile
% imagesc(PNP);axis square
% % nexttile
% % imagesc(PNH);axis square
% % nexttile
% % imagesc(PNR);axis square
% 
% figure;
% tiledlayout(1,2)
% nexttile
% faxis = (-512:511)*Fs/1024;
% ff = fftshift(1:1024);
% ft = abs(fft(UP));
% imagesc(1:1024,faxis,ft(ff,:));axis square
% 
% nexttile
% faxis = (-512:511)*Fs/1024;
% ff = fftshift(1:1024);
% ft = abs(fft(UNP));
% imagesc(1:1024,faxis,ft(ff,:));axis square
% 
% % nexttile
% % faxis = (-512:511)*Fs/1024;
% % ff = fftshift(1:1024);
% % ft = abs(fft(UNH));
% % imagesc(1:1024,faxis,ft(ff,:));axis square
% % 
% % nexttile
% % faxis = (-512:511)*Fs/1024;
% % ff = fftshift(1:1024);
% % ft = abs(fft(UNR));
% % imagesc(1:1024,faxis,ft(ff,:));axis square
%

%% show shifted due to noise
nslv = 0;
Fs = 5e9;
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
PNP = awgn(P,nslv);

[UPN,SPN,VPN] = svd(PNP,'econ');
[UPNf,SPNf,VPNf] = svd(filtfilt(b,a,PNP),'econ');

figure;
tt = tiledlayout(2,5);
for ii=1:5
    nexttile
    plot(taxis,UPN(:,ii),'LineStyle','-','LineWidth',1.5,'Color','Black');
    set(gca,'FontSize',16,'FontName','Times New Roman')
    title(sprintf("U_{%d}",ii),'FontSize',14)
    axis tight;
end

for ii=1:5
    nexttile
    plot(taxis,UPNf(:,ii),'LineStyle','-','LineWidth',1.5,'Color','Black');
    set(gca,'FontSize',16,'FontName','Times New Roman')
    axis tight;
end

xlabel(tt,'Time ($\mu$s)','Interpreter','latex','FontSize',16,'FontName','Times New Roman')
ylabel(tt,'Amplitude','FontSize',16,'FontName','Times New Roman')
%% noise P all level
snrs = linspace(-20,20,20);

% no bpf
dat = zeros(size(snrs));

bpfon = false;
Fs = 5e9;
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
Ploc = P;
if bpfon
    Ploc = filtfilt(b,a,P);
end

[UP,SP,VP] = svd(Ploc,'econ');

for i=1:numel(dat)
    tic
    nslv = snrs(i);
    PNP = awgn(P,nslv);

    if bpfon
        PNP = filtfilt(b,a,PNP);
    end
    
    [UPN,~,~] = svds(PNP,10);
    
    wh = 3;
    normUP = normalize(abs(UP),'range');
    normUPN = normalize(abs(UPN),'range');
    mseval = norm(normUP(:,1:wh)-normUPN(:,1:wh),'fro');
    dat(i) = mseval;
    toc
end
datnbpf = dat;

% with bpf
dat = zeros(size(snrs));

bpfon = true;
Fs = 5e9;
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
Ploc = P;
if bpfon
    Ploc = filtfilt(b,a,P);
end

[UP,SP,VP] = svd(Ploc,'econ');

for i=1:numel(dat)
    tic
    nslv = snrs(i);
    PNP = awgn(P,nslv);

    if bpfon
        PNP = filtfilt(b,a,PNP);
    end
    
    [UPN,~,~] = svds(PNP,10);
    
    wh = 3;
    normUP = normalize(abs(UP),'range');
    normUPN = normalize(abs(UPN),'range');
    mseval = norm(normUP(:,1:wh)-normUPN(:,1:wh),'fro');
    dat(i) = mseval;
    toc
end
datbpf = dat;


%% simul spat dist
T = 1024;
N = 5000;

% H 
l = 128;
nshape = linspace(1,-1,l);
h = [zeros(1,(T-l)/2),nshape,zeros(1,(T-l)/2)];
H = zeros(T,T);
for i=1:T
    H(i,:) = circshift(h,i);
end
H = circshift(H,512,2);

% simulation
pctgs = linspace(1/N,0.25,10);
snrs  = linspace(-20,20,10);

dat = zeros(10,10);

for ii = 1:numel(pctgs)
    for jj = 1:numel(snrs)
        for trial = 1:5
        tic
        fprintf("%d %d\n",ii,jj)
        R = zeros(T,N);
        wh = 3;
        pctg = pctgs(ii);
        lateral_length = round(pctg*N); % consecutive length
        seli = randperm(T,wh);
        selj = randperm(N-round(pctg*N),wh);
        for i=1:wh
            R(seli(i),selj(i):(selj(i)+lateral_length))=1;
        end
        
        P = H*R;
        
        snrval = snrs(jj);% in dB
        % PN = H*awgn(R,snrval);
        PN = awgn(P,snrval);
        %PN = awgn(H,snrval)*R;
        %
        [UP,~,~] = svds(P,wh);
        [UPN,~,~] = svds(PN,wh);
        
        normUP = normalize(abs(UP),'range');
        normUPN = normalize(abs(UPN),'range');
        
        mseval = norm(normUP(:,1:wh)-normUPN(:,1:wh),'fro');
        dat(ii,jj) = dat(ii,jj)+mseval;
        toc
        end
    end
end
%%
pctgs = linspace(1/N,0.25,10);
snrs  = linspace(-20,20,10);
PNP = awgn(P,5);
[UPN,SPN,VPN] = svd(PNP,'econ');
[UPNf,SPNf,VPNf] = svd(filtfilt(b,a,PNP),'econ');

figure;
tt = tiledlayout(1,3);

nexttile
faxis = (-512:511)*Fs/1024;
ff = fftshift(1:1024);
ft = mag2db(abs(fft(UPN)));
imagesc(1:1024,faxis/1e9,ft(ff,:));hold off
hcb=colorbar;
hcb.Title.String = "(dB)";
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("$|\mathcal{F} (U_P)|$",'Interpreter', 'LaTeX')
xlabel("# Number of singular values",'FontSize',16,'FontName','Times New Roman')
ylabel("Frequency (GHz)",'FontSize',16,'FontName','Times New Roman')

nexttile
faxis = (-512:511)*Fs/1024;
ff = fftshift(1:1024);
ft = mag2db(abs(fft(UPNf)));
imagesc(1:1024,faxis/1e9,ft(ff,:));hold off
hcb=colorbar;
hcb.Title.String = "(dB)";
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("$|\mathcal{F} (U_{P_{BPF}})|$",'Interpreter', 'LaTeX')

xlabel("# Number of singular values",'FontSize',16,'FontName','Times New Roman')
ylabel("Frequency (GHz)",'FontSize',16,'FontName','Times New Roman')

nexttile
imagesc(snrs,pctgs,dat/5);
xlabel("SNR (dB)")
ylabel("$k$", Interpreter='latex')
set(gca,'FontName','Times','FontSize',16)
axis square;colormap jet;
hcb=colorbar;
hcb.Title.String = "m.s.e";
%% make sure the awgn works correctly
ns = awgn(zeros(1,T),0);
figure;
plot(mag2db(abs(fft(h))));hold on;
plot(mag2db(abs(fft(ns))));hold off;
axis tight