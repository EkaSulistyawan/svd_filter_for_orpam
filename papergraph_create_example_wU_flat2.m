%% load data 
close all
load("../../SharingPoint/data/cellular/compiled/04092023_USAF.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:100,1:100);
Fs = 5e9;
%% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);


%% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,100*100);
[u,s,v] = svd(dat2d,"econ");

%% show example of wU
sig = u;
minsig = repmat(min(sig),[size(u,2) 1]);
rangesig = repmat(range(sig),[size(u,2) 1]);
sig = (sig - minsig) ./ rangesig;
%normalized = (sig*2 - 1.0);
normalized = (sig - 0.5)*2;
energy = (normalized.^2);

wU = -sum(energy);


%%
timeaxis = (1:1024)*1e9/5e9;
figure('Position',[592,194,700,746])
subplot(321);
plot(timeaxis,sig(:,1));
title("Signal");
ylabel("Intensity (a.u.)","FontWeight","bold");
set(gca,'XTick',[],'FontSize',15,'FontName','Times New Roman');

subplot(322);
plot(timeaxis,sig(:,9));
title("Noise");
set(gca,'XTick',[],'YTick',[],'FontSize',15,'FontName','Times New Roman');

subplot(323);
plot(timeaxis,normalized(:,1));
ylabel("Normalized [-1,1]","FontWeight","bold")
set(gca,'XTick',[],'FontSize',15,'FontName','Times New Roman');

subplot(324);
plot(timeaxis,normalized(:,9));
set(gca,'XTick',[],'YTick',[],'FontSize',15,'FontName','Times New Roman');

subplot(325);
a = area(timeaxis,energy(:,1));
a(1).FaceColor = [0.8 0.8 0.8];
a(1).FaceAlpha = 0.5;
a(1).EdgeColor = [0 0.4470 0.7410];
%plot(timeaxis,energy(:,1));
ylabel("Squared (0,1]","FontWeight","bold")
xlabel("Time (ns)")
set(gca,'FontSize',15,'FontName','Times New Roman');

subplot(326);
a = area(timeaxis,energy(:,9));
a(1).FaceColor = [0.8 0.8 0.8];
a(1).FaceAlpha = 0.5;
a(1).EdgeColor = [0 0.4470 0.7410]; % Matlab's default
%plot(timeaxis,energy(:,2));
xlabel("Time (ns)")
set(gca,'YTick',[],'FontSize',15,'FontName','Times New Roman');

exportgraphics(gcf,'./exported_images/fig5_wU_flat2.emf')