%% load data 
close all
load("../../SharingPoint/data/cellular/compiled/17122022_USAF.mat")
% load("../data/cellular/compiled/17122022_USAF.mat");
space3D = space3D(:,1:200,1:200);
Fs = 5e9;
%% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);


%% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,200*200);
[u,s,v] = svd(dat2d,"econ");

%% create the cmode anddata2D

figure('Position',[672,427,1106,420])
ax(1) = subplot(121);
cmode_func = @(x)(squeeze(max(abs(hilbert(x)))));
imagesc(cmode_func(space3Dbpf)');
hcb = colorbar;
ylabel(hcb,'Intensity (a.u.)',"FontName",'Times','FontSize',18)
title("(A)")
set(gca,"FontName",'Times','FontSize',18)
line([156 176],[188 188],'LineWidth',5,'Color','white')
axis off

% ax = gca;
% exportgraphics(ax,'./exported_images/fig1_cmode.jpg','Resolution',900)

ax(2) = subplot(122);
taxis = (1:1024)*1e9/Fs;
xaxis = 1:1024;
imagesc(xaxis,taxis,dat2d);
xlabel("Number of A-lines (#)")
ylabel("Time (ns)")
title("(B)")
hcb = colorbar;
ylabel(hcb,'Intensity (a.u.)',"FontName",'Times','FontSize',18)
set(gca,"FontName",'Times','FontSize',18)

colormap(ax(1),hot)
colormap(ax(2),gray)
ax = gcf;
exportgraphics(ax,'./exported_images/fig3_cmode_2d.emf')

%% create the U, S, V
figure('Position',[13 422 1534 420])

subplot(131)
taxis = (1:1024)*1e9/Fs;
xaxis = 1:1024;
imagesc(xaxis,taxis,u);
xlabel("Number of A-lines (#)")
ylabel("Time (ns)")
title("(C)")
%hcb = colorbar;
%ylabel(hcb,'Intensity (a.u.)',"FontName",'Times','FontSize',18)
set(gca,"FontName",'Times','FontSize',18)
colormap gray

subplot(132)
taxis = (1:1024)*1e9/Fs;
xaxis = 1:1024;
plot(xaxis,diag(s));
xlabel("Number of A-lines (#)")
ylabel("Magnitude (a.u.)")
title("(D)")
set(gca,"FontName",'Times','FontSize',18)
colormap gray

subplot(133)
imagesc(v);
xlabel("Number of A-lines (#)")
ylabel("Spatial Vectors (a.u.)")
title("(E)")
%hcb = colorbar;
%ylabel(hcb,'Intensity (a.u.)',"FontName",'Times','FontSize',18)
set(gca,"FontName",'Times','FontSize',18)
colormap gray

ax = gcf;
exportgraphics(ax,'./exported_images/fig3_usv.emf')
%% Individual; get the first C-mode, B-mode and the signal
sel = [1:10];
% see the A line
subplot(3,1,1);plot(sum(u(:,sel),2));
xlim([1, 1024])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,200,200);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,1,2);imagesc(cmode);axis off;colormap hot

bmode = squeeze(hb(:,99,:));
subplot(3,1,3);imagesc(bmode);axis off;colormap hot

%% Collection; get the first C-mode, B-mode and the signal
numSel = 6;
axesCount = 1;

figure('Position',[1,41,1920,963])
for p=1:numSel
sel = p;
% see the A line
maxval = max(u(:,1:numSel),[],"all");
minval = min(u(:,1:numSel),[],"all");
taxis = (1:1024)*1e9/Fs;
subplot(3,numSel,axesCount);plot(taxis,u(:,sel));
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([minval, maxval])
xlim([1, max(taxis)])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,200,200);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,numSel,axesCount+numSel);
imagesc(cmode');axis off;colormap hot
if(p == 1)
    line([1 200],[32 32],'LineWidth',2,'Color','g','LineStyle','--')
end

bmode = squeeze(hb(:,:,32));
xaxis = (1:200)*200 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
subplot(3,numSel,axesCount+2*numSel);imagesc(xaxis,yaxis,bmode);colormap hot;
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])

axesCount = axesCount+1;
end
ax = gcf;
exportgraphics(ax,'./exported_images/fig4_exampleUSAF.emf')