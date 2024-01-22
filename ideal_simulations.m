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
rule3 = randperm(N,1);
rule3 = rule3:(rule3+0.04*N);
R(200,rule3) = 1;

% rule 4
rule4 = randperm(N,1);
rule4 = rule4 :(rule4+0.08*N);
R(800,rule4) = 1;

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
leg = legend('Orientation', 'Horizontal','FontSize',14);

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
    plot(conv(abs(eq5b(:,sel)),mav),'Color',[0 0 0 0.2],'DisplayName','$R^T V_H$','LineWidth',1.5);hold on
    plot(conv(abs(VR(:,sel)),mav),'Color',[1 0 1 0.2],'DisplayName','$V_R$','LineWidth',1.5);hold on
    plot(conv(abs(VP(:,sel)),mav),'Color',[1 0 0 0.2],'DisplayName','$V_P$','LineWidth',1.5);hold off
    set(gca,'FontSize',16,'FontName','Times New Roman')
    axis square tight
    ylim([0 1])
end
legend('Interpreter', 'LaTeX')
leg = legend('Orientation', 'Horizontal','FontSize',14);

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