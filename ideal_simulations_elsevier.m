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
axis square tight
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
% ylabel("Time ($\mu$s)",'Interpreter','latex')
axis image
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it H")


% construct R
R = zeros(T,N);


% rule 2
% rule2 = 1:N;
% rule3 = [1:1000,2000:4000,4500:5000];
% R(512,rule3) = 1;

% rule 3
% rule3 = randperm(N,1);
% rule3 = rule3:(rule3+0.04*N);
% R(200,rule3) = 1;

% rule 4
% rule4 = randperm(N,1);
% rule4 = rule4 :(rule4+0.08*N);
% R(800,rule4) = 1;

% rule 5
rule1 = 1000:3000;
R(200,rule1) = 1;

% rule 6
rule2 = 4200:4500;
R(800,rule2) = 1;

rule3 = 100:150;
R(500,rule3) = 1;

Rbefrand = R;
kpctg = 1/T;
rnd = randperm(T*N,kpctg*T*N);
% R(rnd) = 1;

nexttile
% for visualization purpose
Rvis = Rbefrand;
for i=1:20
    Rvis = Rvis + circshift(R,i,1);
end
Rvis(rnd) = 1;
Rvis(Rvis > 1) = 1;
imagesc(xaxis,taxis,Rvis);
xlabel(sprintf("Lateral \nPosition (a.u.)"))
%ylabel("Time ($\mu$s)",'Interpreter','latex')
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it R")


% construct P
P = H*R;

nexttile
imagesc(xaxis,taxis,P);
xlabel(sprintf("Lateral \nPosition (a.u.)"))
ylabel("Time ($\mu$s)",'Interpreter','latex')
axis square
set(gca,'FontSize',16,'FontName','Times New Roman')
title("\it P")

%% compute SVD
[UP,SP,VP] = svd(P);
[UH,SH,VH] = svd(H);
[UR,SR,VR] = svd(R);






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
        wh = 2;
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

%% overlapping artifact
close all
% exectue the first section first to get P, H and R
naxis = linspace(0,1,1024);
xi = zeros(T,N);
pos =1:2:5000; 
artf = sin(2*pi*2*naxis);
xi(:,pos) = 0.1*repmat(artf',1,numel(pos));

fntsz = 20;

% rotate
sine(:,1) = taxis;
sine(:,2) = artf;

theta = pi/2;
RotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
sine2 = [sine(:,1),sine(:,2)] * RotationMatrix;

Pxi = P + xi;
Pxi = awgn(Pxi,0);


[up,sp,vp] = svd(P,'econ');
[upxi,spxi,vpxi] = svd(Pxi,'econ');

figure
tiledlayout(1,4)

tt1=nexttile;
imagesc(xaxis,taxis,P);
colormap(tt1,'gray')
xlabel(sprintf("Lateral \nPosition (a.u.)"))
ylabel("Time ($\mu$s)",'Interpreter','latex')
axis square tight
set(gca,'FontSize',fntsz,'FontName','Times New Roman')
title("$\it{P}$",'Interpreter','latex')

tt2=nexttile;
imagesc(xaxis,taxis,Pxi); hold on
plot(sine2(:,1)*100+200,sine2(:,2)+taxis(end),'Color','green'); hold off
colormap(tt2,'gray')
xlabel(sprintf("Lateral \nPosition (a.u.)"))
ylabel("Time ($\mu$s)",'Interpreter','latex')
axis square tight
set(gca,'FontSize',fntsz,'FontName','Times New Roman')
title("$\it{P}+\xi+\eta$",'Interpreter','latex')


tt3 = nexttile;
alp = 0.3;
colorset=[[1 0 0 alp];[0 0 0 alp];[0 0 1 alp]];
for i=1:3
    %plot(up(:,i),'LineStyle','-','LineWidth',1.5,'Color','Black','DisplayName','H U_R');hold on
    plot(tt3,taxis,upxi(:,i),'LineStyle','-','LineWidth',1.5,'Color',colorset(i,:),'DisplayName',sprintf("#%d",i));hold on
    axis square tight
    set(gca,'FontSize',fntsz,'FontName','Times New Roman')
end
xlabel("Time ($\mu$s)",'Interpreter','latex')
legend('Orientation','horizontal','location','south')
title('$U_{\it{P}+\xi+\eta}$','Interpreter','Latex')

tt4 = nexttile;
% use the 3 consecutive position
load("experiment28022024.mat",'dat');
pctgs = linspace(1/N,0.25,10);
snrs  = linspace(-20,20,10);
imagesc(snrs,ceil(pctgs*100),dat/5);
colormap(tt4,'jet')
xlabel("SNR (dB)")
ylabel("$k\%$", Interpreter='latex')
set(gca,'FontName','Times','FontSize',fntsz)
axis square;
hcb=colorbar;
hcb.Title.String = "m.s.e";
title('MSE of $U_{\it{P}+\eta}$','Interpreter','Latex')


annotation('textbox',[0.1 0.1 0.1 0.1],'String','$k=6\%$','EdgeColor','none','FontSize',fntsz,'FontName','Times','Color','green','Interpreter','latex')
annotation('textbox',[0.1 0.1 0.1 0.1],'String','$\xi$','EdgeColor','none','FontSize',fntsz,'FontName','Times','Color','green','Interpreter','latex')

%% make sure the awgn works correctly
ns = awgn(zeros(1,T),0);
figure;
plot(mag2db(abs(fft(h))));hold on;
plot(mag2db(abs(fft(ns))));hold off;
axis tight