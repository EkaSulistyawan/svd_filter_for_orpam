%% flat
load("../../SharingPoint/data/cellular/compiled/04092023_USAF.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:100,1:100);
Fs = 5e9;
% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);

% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,100*100);
[u,s,v] = svd(dat2d,"econ");

% show example of wU
sig = u;
minsig = repmat(min(sig),[size(u,2) 1]);
rangesig = repmat(range(sig),[size(u,2) 1]);
sig = (sig - minsig) ./ rangesig;
sig = (sig - 0.5)*2;
energy = -sum(sig.^2);
energy = (energy - min(energy))./ range(energy);
%wU = energy' .* diag(s);
energy_flat = energy;

imsize = size(space3D,2);
v3d = reshape(v,imsize,imsize,1024);
objectness = get_dctweight(v3d);
objectness = (objectness - min(objectness)) ./ range(objectness);
%wV = objectness' .* diag(s);
objectness_flat = objectness;


% gamma=0  ; w0 = (gamma*wU) + (1-gamma)*wV;
% gamma=0.5; whalf = (gamma*wU) + (1-gamma)*wV;
% gamma=1  ; w1 = (gamma*wU) + (1-gamma)*wV;

%%
figure
sel = [1:20];
plot(sel,energy(sel),'DisplayName','w_U','Marker','+');hold on
%plot(sel,whalf(sel),'DisplayName','w\Sigma (\gamma=0.5)','Marker','o');hold on
plot(sel,objectness(sel),'DisplayName','w_V','Marker','x');hold off
xlim([1 20])
legend('Location','southwest')

set(gca,'XTick',[1:20],'XTickLabel',[1:20],'FontSize',15,'FontName','Times New Roman');


%% round
load("../../SharingPoint/data/cellular/compiled/09102022_sphere.mat")
space3D = space3D(:,1:200,1:200);
Fs = 5e9;
% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);

% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,200*200);
[u,s,v] = svd(dat2d,"econ");

% show example of wU
sig = u;
minsig = repmat(min(sig),[size(u,2) 1]);
rangesig = repmat(range(sig),[size(u,2) 1]);
sig = (sig - minsig) ./ rangesig;
sig = (sig - 0.5)*2;
energy = -sum(sig.^2);
energy = (energy - min(energy))./ range(energy);
%wU = energy' .* diag(s);
energy_round = energy;


imsize = size(space3D,2);
v3d = reshape(v,imsize,imsize,1024);
objectness = get_dctweight(v3d);
objectness = (objectness - min(objectness)) ./ range(objectness);
%wV = objectness' .* diag(s);
objectness_round = objectness;


% gamma=0  ; w0 = (gamma*wU) + (1-gamma)*wV;
% gamma=0.5; whalf = (gamma*wU) + (1-gamma)*wV;
% gamma=1  ; w1 = (gamma*wU) + (1-gamma)*wV;

%%
figure('Position',[463,412,1307,420])
tiledlayout(1,2)
nexttile
sel = [1:20];
plot(sel,energy_flat(sel),'DisplayName','w_U','Marker','+');hold on
%plot(sel,whalf(sel),'DisplayName','w\Sigma (\gamma=0.5)','Marker','o');hold on
plot(sel,objectness_flat(sel),'DisplayName','w_V','Marker','x');hold off
title("Weight for USAF-1951")
ylabel("Weight Amplitude [0 1]")
xlabel("Singular Vectors")
xlim([1 20])
legend('Location','southwest')
set(gca,'XTick',[1:2:20],'XTickLabel',[1:2:20],'FontSize',15,'FontName','Times New Roman');

nexttile
plot(sel,energy_round(sel),'DisplayName','w_U','Marker','+');hold on
%plot(sel,whalf(sel),'DisplayName','w\Sigma (\gamma=0.5)','Marker','o');hold on
plot(sel,objectness_round(sel),'DisplayName','w_V','Marker','x');hold off
title("Weight for Red-stained Spherical Phantom")
ylabel("Weight Amplitude [0 1]")
xlabel("Singular Vectors")
xlim([1 20])
legend('Location','southwest')
set(gca,'XTick',[1:2:20],'XTickLabel',[1:2:20],'FontSize',15,'FontName','Times New Roman');

exportgraphics(gcf,'./exported_images/weight.emf')