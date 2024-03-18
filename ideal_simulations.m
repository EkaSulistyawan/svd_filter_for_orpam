clear all 
close all


% Simulation specs
N = 5000;
T = 1024;
l = 128;   % PA length


tiledlayout(1,4);
% Construct h
nshape = linspace(0.5,-0.5,l);
h = [zeros(1,(T-l)/2),nshape,zeros(1,(T-l)/2)];

% fs = 5e9;
% [b,a] = butter(4,[10e6 100e6] / (fs/2),'bandpass');
% h = filtfilt(b,a,h);

nexttile
plot(1:T,h,'LineWidth',2,'Color','black');
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("{\it h(\it{t})}")

% Construct H
H = zeros(T,T);
for i=1:T
    H(i,:) = circshift(h,i);
end
H = circshift(H,512,2);

nexttile
imagesc(H);
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
imagesc(R);
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it R")

% construct P
P = H*R;

nexttile
imagesc(P);
axis square
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

%% Proof 2 relation UP and UH

figure;
tiledlayout(1,2)
nexttile
imagesc(abs(fft(UP)));hold on
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("$|\mathcal{F} (U_P)|$",'Interpreter', 'LaTeX')

nexttile
imagesc(abs(fft(UH)));hold off
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("$|\mathcal{F} (U_H)|$",'Interpreter', 'LaTeX')


%%
fewranks = 4;
eq5b = R'*VH;

mav = ones(1,10);

figure;
tiledlayout(1,fewranks)
for sel =1:fewranks
    nexttile
    %plot(conv(abs(eq5b(:,sel)),mav),'Color',[0 0 0 0.2],'DisplayName','$R^T V_H$','LineWidth',1.5);hold on
    plot(conv(abs(VR(:,sel)),mav),'Color',[0 0 0 0.8],'DisplayName','$V_R$','LineWidth',1.5);hold on
    plot(conv(abs(VP(:,sel)),mav),'Color',[1 0 0 0.8],'DisplayName','$V_P$','LineWidth',1.5);hold off
    set(gca,'FontSize',16,'FontName','Times New Roman')
    axis square tight
    ylim([0 1])
end
% legend('Interpreter', 'LaTeX')
% leg = legend('Orientation', 'Horizontal','FontSize',14);

%% sigma
figure;
semilogy(diag(SP),'DisplayName','$\Sigma_P$','LineWidth',1.5);hold on;
semilogy(diag(SH),'DisplayName','$\Sigma_H$','LineWidth',1.5);hold on;
semilogy(diag(SR),'DisplayName','$\Sigma_R$','LineWidth',1.5);hold off;
legend('Interpreter','latex','Location','southeast')
set(gca,'FontSize',16,'FontName','Times New Roman')
axis tight square

axes('position',[.3 .4 .3 .3])
box on

tmp = diag(SP);semilogy(tmp(1:20),'DisplayName','$\Sigma_P$','LineWidth',1.2);hold on;
tmp = diag(SH);semilogy(tmp(1:20),'DisplayName','$\Sigma_H$','LineWidth',1.2);hold on;
tmp = diag(SR);semilogy(tmp(1:20),'DisplayName','$\Sigma_R$','LineWidth',1.2);hold off;
set(gca,'FontSize',12,'FontName','Times New Roman')

%% sigma 2
figure;
semilogy(diag(SP),'DisplayName','$\Sigma_P$','LineWidth',1.5);hold on;
semilogy(diag(SH).*diag(SR),'DisplayName','$\Sigma_H \Sigma_R$','LineWidth',1.5);hold off;
legend('Interpreter','latex','Location','southeast')
set(gca,'FontSize',16,'FontName','Times New Roman')
axis tight square

axes('position',[.3 .4 .3 .3])
box on

tmp = diag(SP);semilogy(tmp(1:20),'DisplayName','$\Sigma_P$','LineWidth',1.2);hold on;
tmp = diag(SR).*diag(SH);semilogy(tmp(1:20),'DisplayName','$\Sigma_H \Sigma_R$','LineWidth',1.2);hold off;
set(gca,'FontSize',12,'FontName','Times New Roman')

%% recover R from H
% get H 
h2 = UP(:,1);
H2 = zeros(T,T);
for i=1:T
    H2(i,:) = circshift(h2,i);
end
H2 = circshift(H2,512,2);

% based on Eq 5a

% based on Eq 5b
VR2 = VP(:,1:3);

%% UP filt
[u ,s ,v] = svd(P,'econ');
[uf,sf,vf] = svd(fft(P),'econ');
urec = real(ifft(uf));
kk = real(ifft(uf*sf*vf'));

%% effect of noise to P
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

%% add noise to P
noised = awgn(P,-10);
[UN,SN,VN] = svd(noised);
figure;plot(UN(:,1:1))
figure;imagesc(noised)


%% add noise to H
PN = awgn(H,20)*R;
[UN,SN,VN] = svd(PN);
figure;plot(UN(:,1:3))
figure;
tiledlayout(1,2)
nexttile
imagesc(P);axis square
nexttile
imagesc(PN);axis square

figure;
tiledlayout(1,2)
nexttile
imagesc(abs(fft(UP)));axis square

nexttile
imagesc(abs(fft(UN)));axis square

%% add noise to R
PN = H*awgn(R,10);
[UN,SN,VN] = svd(PN);
figure;plot(UN(:,1:3))
figure;
tiledlayout(1,2)
nexttile
imagesc(P);axis square
nexttile
imagesc(PN);axis square

figure;
tiledlayout(1,2)
nexttile
imagesc(abs(fft(UP)));axis square

nexttile
imagesc(abs(fft(UN)));axis square

%%
figure;
tiledlayout(1,2)

%% study of R
% close all
% since R holds the strongest influence towards the SVD, we may consider
% how R should be constructed when noise or inhomogenousity existed.

% at its extreme, the signal is at the same level as the inhomogenousity.
% Ofc in this case, we can only rely on spatial coherence.

pctgs = linspace(0.05,0.1,10);
nsadds = linspace(0,15,10);

dat = zeros(numel(pctgs),numel(nsadds));

colnum = 1000;

for i = 1:numel(pctgs)
    for j = 1:numel(nsadds)
        pctg = pctgs(i);
        nsadd = nsadds(j);
        % new R
        sglv = 1;
        R2 = zeros(1024,colnum);
        % pctg = 0.05;
        % R2(512,1:round(pctg*colnum)) = sglv;
        R2(512,randperm(colnum,round(pctg*colnum))) = sglv;
        % R2(10,randperm(5000,100)) = 1;
        % R2(100,randperm(5000,round(pctg*5000))) = sglv;
        
        
        % old R
        %R2 = R;
        snrlv = -10*log(sglv)+nsadd;
        Rn = awgn(R2,snrlv);
        % fprintf("Rank Noised %d\n",rank(Rn))
        % [UR,SR,VR] = svd(R2);
        [URn,SRn,VRn] = svd(Rn,'econ');
        
        Un = URn(:,1);
        Un(512) = 0;
        varnoise = std(Un);
        significantEnergy = URn(512,1) - varnoise;

        dat(i,j) = significantEnergy;
    end
end

figure;
imagesc(nsadds,pctgs,abs(dat));
axis square;colormap jet;
colorbar
set(gca,'FontSize',20)
clim([0 1])
xlabel("SNR")
ylabel("% Number of Absorber")

% figure;
%semilogy(diag(SR),'DisplayName','R');hold on
% semilogy(diag(SRn),'DisplayName',sprintf('R noisy. Signal level %d dB',snrlv));hold on
% legend
% axis tight

% % figure;
% plot(URn(:,1)); hold on
% axis tight
% 
% figure;imagesc(Rn*Rn');axis square;colormap gray
% figure;imagesc(Rn'*Rn);axis square;colormap gray

% figure;plot(diag(Rn*Rn'));axis tight;colormap gray
% figure;plot(diag(Rn'*Rn));axis tight;colormap gray
% what are noise level? 

