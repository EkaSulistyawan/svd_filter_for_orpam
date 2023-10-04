%% load data 
close all
clear all
%load("../data/cellular/compiled/11102022_sphere2.mat");
load("../../SharingPoint/data/cellular/compiled/09102022_sphere.mat")
imsize = size(space3D,2)-1;
space3D = space3D(:,1:imsize,1:imsize);
Fs = 5e9;

%% BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space3Dbpf = filtfilt(b,a,space3D);

%% SVD on bandpassed filter data
dat2d = reshape(space3Dbpf,1024,imsize*imsize);
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
exportgraphics(ax,'./exported_images/fig5_cmode_2d.emf')

%% Individual; get the first C-mode, B-mode and the signal
sel = [1:10];
% see the A line
subplot(3,1,1);plot(sum(u(:,sel),2));
xlim([1, 1024])

% reconstruct the data
srec = zeros(size(s));
srec(sel,sel) = s(sel,sel);
rec = u*srec*v';
rec3d= reshape(rec,1024,imsize,imsize);

% get the cmode image
hb = abs(hilbert(rec3d));
cmode = squeeze(max(hb));
subplot(3,1,2);imagesc(cmode);axis off;colormap hot

bmode = squeeze(hb(:,99,:));
subplot(3,1,3);imagesc(bmode);axis off;colormap hot



%% Collection; get the first C-mode, B-mode and the signal
% numSel = 6;
% axesCount = 1;
% 
% figure('Position',[1,41,1920,963])
% for p=1:numSel
% sel = p;
% % see the A line
% maxval = max(u(:,1:numSel),[],"all");
% minval = min(u(:,1:numSel),[],"all");
% taxis = (1:1024)*1e9/Fs;
% subplot(3,numSel,axesCount);plot(taxis,u(:,sel));
% xlabel("Time (ns)")
% ylabel("Intensity (a.u.)")
% ylim([minval, maxval])
% xlim([1, max(taxis)])
% 
% % reconstruct the data
% srec = zeros(size(s));
% srec(sel,sel) = s(sel,sel);
% rec = u*srec*v';
% rec3d= reshape(rec,1024,200,200);
% 
% % get the cmode image
% hb = abs(hilbert(rec3d));
% cmode = squeeze(max(hb));
% subplot(3,numSel,axesCount+numSel);
% imagesc(cmode');axis off;colormap hot
% if(p == 1)
%     line([1 200],[159 159],'LineWidth',2,'Color','g','LineStyle','--')
% end
% 
% bmode = squeeze(hb(:,:,159));
% xaxis = (1:200)*200 / 1000; 
% yaxis = (1:1024)*1e6*1480/Fs;
% subplot(3,numSel,axesCount+2*numSel);imagesc(xaxis,yaxis,bmode);colormap hot;
% xlabel(['Length (' char(181) 'm)'])
% ylabel(['Depth (' char(181) 'm)'])
% 
% axesCount = axesCount+1;
% end
% sgtitle("From 1 to 6",'FontName','Times')
% ax = gcf;
% exportgraphics(ax,'./exported_images/fig5_exampleSphereA.emf')

%% Collection; get the first C-mode, B-mode and the signal
% numSel = 6;
% axesCount = 1;
% 
% figure('Position',[1,41,1920,963])
% for p=1:numSel
% sel = p + 6;
% % see the A line
% maxval = max(u(:,1:numSel),[],"all");
% minval = min(u(:,1:numSel),[],"all");
% taxis = (1:1024)*1e9/Fs;
% subplot(3,numSel,axesCount);plot(taxis,u(:,sel));
% xlabel("Time (ns)")
% ylabel("Intensity (a.u.)")
% ylim([minval, maxval])
% xlim([1, max(taxis)])
% 
% % reconstruct the data
% srec = zeros(size(s));
% srec(sel,sel) = s(sel,sel);
% rec = u*srec*v';
% rec3d= reshape(rec,1024,200,200);
% 
% % get the cmode image
% hb = abs(hilbert(rec3d));
% cmode = squeeze(max(hb));
% subplot(3,numSel,axesCount+numSel);
% imagesc(cmode');axis off;colormap hot
% if(p == 1)
%     line([1 200],[159 159],'LineWidth',2,'Color','g','LineStyle','--')
% end
% 
% bmode = squeeze(hb(:,:,159));
% xaxis = (1:200)*200 / 1000; 
% yaxis = (1:1024)*1e6*1480/Fs;
% subplot(3,numSel,axesCount+2*numSel);imagesc(xaxis,yaxis,bmode);colormap hot;
% xlabel(['Length (' char(181) 'm)'])
% ylabel(['Depth (' char(181) 'm)'])
% 
% axesCount = axesCount+1;
% end
% %sgtitle("From 7 to 12",'FontName','Times')
% ax = gcf;
% exportgraphics(ax,'./exported_images/fig5_exampleSphereB.emf')

%% Collection; get the first C-mode, B-mode and the signal
% numSel = 6;
% axesCount = 1;
% 
% figure;
% for p=1:numSel
% sel = p + numSel;
% % see the A line
% maxval = max(u(:,1:numSel),[],"all");
% minval = min(u(:,1:numSel),[],"all");
% subplot(3,numSel,axesCount);plot(u(:,sel));
% ylim([minval, maxval])
% xlim([1, 1024])
% 
% % reconstruct the data
% srec = zeros(size(s));
% srec(sel,sel) = s(sel,sel);
% rec = u*srec*v';
% rec3d= reshape(rec,1024,imsize,imsize);
% 
% % get the cmode image
% hb = abs(hilbert(rec3d));
% cmode = squeeze(max(hb));
% subplot(3,numSel,axesCount+numSel);imagesc(cmode);axis off;colormap hot
% 
% bmode = squeeze(hb(:,36,:));
% subplot(3,numSel,axesCount+2*numSel);imagesc(bmode);axis off;colormap hot
% 
% axesCount = axesCount+1;
% end

%% Collection; get the first C-mode, B-mode and the signal
% numSel = 6;
% axesCount = 1;
% 
% figure;
% for p=1:numSel
% sel = p + 2*numSel;
% % see the A line
% maxval = max(u(:,1:numSel),[],"all");
% minval = min(u(:,1:numSel),[],"all");
% subplot(3,numSel,axesCount);plot(u(:,sel));
% ylim([minval, maxval])
% xlim([1, 1024])
% 
% % reconstruct the data
% srec = zeros(size(s));
% srec(sel,sel) = s(sel,sel);
% rec = u*srec*v';
% rec3d= reshape(rec,1024,imsize,imsize);
% 
% % get the cmode image
% hb = abs(hilbert(rec3d));
% cmode = squeeze(max(hb));
% subplot(3,numSel,axesCount+numSel);imagesc(cmode);axis off;colormap hot
% 
% bmode = squeeze(hb(:,36,:));
% subplot(3,numSel,axesCount+2*numSel);imagesc(bmode);axis off;colormap hot
% 
% axesCount = axesCount+1;
% end

%% try tiledlayout
numSel = 6;
selectedSel = [1,2,3,10,11,12];
%selectedSel = [1:6];
totalCol = 6;
axesCount = 1;

% find the highest maxval
% find the minimum minval
maxval = max(u(:,1:numSel),[],"all");
minval = min(u(:,1:numSel),[],"all");
figure('Position',[1955,167,1887,849])
mytile = tiledlayout(3,numSel);
for p=1:(numSel)
    sel = selectedSel(p);
    
    % see the A line
    nexttile
    taxis = (1:1024)*1e9/Fs;
    plot(taxis,u(:,sel));
    xlabel("Time (ns)")
    ylabel("Intensity (a.u.)")
    ylim([minval, maxval])
    xlim([1, max(taxis)])
    set(gca,"FontName","Times New Roman","FontSize",14)
end

%
for p=1:(numSel)
    sel = selectedSel(p);
    % reconstruct the data
    srec = zeros(size(s));
    srec(sel,sel) = s(sel,sel);
    rec = u*srec*v';
    rec3d= reshape(rec,1024,200,200);
    
    % get the cmode image
    nexttile(p + totalCol)
    hb = abs(hilbert(rec3d));
    cmode = squeeze(max(hb));
    imagesc(cmode');axis off;colormap hot;
    if(p == 1)
        line([1 200],[159 159],'LineWidth',2,'Color','g','LineStyle','--')
        line([156 176],[100 100],'LineWidth',2,'Color','white');
        
    end
end


for p=1:numSel
    sel = selectedSel(p);
    % reconstruct the data
    srec = zeros(size(s));
    srec(sel,sel) = s(sel,sel);
    rec = u*srec*v';
    rec3d= reshape(rec,1024,200,200);
    hb = abs(hilbert(rec3d));

    % get the bmode
    nexttile(p + 2*totalCol)
    bmode = squeeze(hb(:,:,159));
    xaxis = (1:200)*400 / 1000; 
    yaxis = (1:1024)*1e6*1480/Fs;
    imagesc(xaxis,yaxis,bmode);colormap hot;
    xlabel(['Length (' char(181) 'm)'])
    ylabel(['Depth (' char(181) 'm)'])
    set(gca,"FontName","Times New Roman","FontSize",14)

end

annotation('arrow',[0.12 0.1],[0.6 0.6],'Color','g');
annotation('arrow',[0.42 0.4],[0.6 0.6],'Color','g');

ax = gcf;
exportgraphics(ax,'./exported_images/fig4_exampleSphere.emf')

%% try tiledlayout
numSel = 6;
%selectedSel = [1,2,3,10,11,12];
selectedSel = [7:12];
totalCol = 6;
axesCount = 1;

% find the highest maxval
% find the minimum minval
maxval = max(u(:,1:numSel),[],"all");
minval = min(u(:,1:numSel),[],"all");
figure('Position',[1955,167,1887,849])
mytile = tiledlayout(3,numSel);
for p=1:(numSel)
    sel = selectedSel(p);
    
    % see the A line
    nexttile
    taxis = (1:1024)*1e9/Fs;
    plot(taxis,u(:,sel));
    xlabel("Time (ns)")
    ylabel("Intensity (a.u.)")
    ylim([minval, maxval])
    xlim([1, max(taxis)])
    set(gca,"FontName","Times New Roman","FontSize",14)
end

%
for p=1:(numSel)
    sel = selectedSel(p);
    % reconstruct the data
    srec = zeros(size(s));
    srec(sel,sel) = s(sel,sel);
    rec = u*srec*v';
    rec3d= reshape(rec,1024,200,200);
    
    % get the cmode image
    nexttile(p + totalCol)
    hb = abs(hilbert(rec3d));
    cmode = squeeze(max(hb));
    imagesc(cmode');axis off;colormap hot;
    if(p == 1)
        line([1 200],[159 159],'LineWidth',2,'Color','g','LineStyle','--')
        line([156 176],[100 100],'LineWidth',2,'Color','white');
        
    end
end


for p=1:numSel
    sel = selectedSel(p);
    % reconstruct the data
    srec = zeros(size(s));
    srec(sel,sel) = s(sel,sel);
    rec = u*srec*v';
    rec3d= reshape(rec,1024,200,200);
    hb = abs(hilbert(rec3d));

    % get the bmode
    nexttile(p + 2*totalCol)
    bmode = squeeze(hb(:,:,159));
    xaxis = (1:200)*400 / 1000; 
    yaxis = (1:1024)*1e6*1480/Fs;
    imagesc(xaxis,yaxis,bmode);colormap hot;
    xlabel(['Length (' char(181) 'm)'])
    ylabel(['Depth (' char(181) 'm)'])
    set(gca,"FontName","Times New Roman","FontSize",14)

end

%annotation('arrow',[0.12 0.1],[0.6 0.6],'Color','g');
%annotation('arrow',[0.42 0.4],[0.6 0.6],'Color','g');

ax = gcf;
exportgraphics(ax,'./exported_images/exampleSphere_7to12.emf')

%%
hf = figure('Units','normalized'); 
colormap hot
cbh = colorbar;
cbh.Location = "south";
cbh.Label.String="Normalized Intensity [0, 1]";
set(gca,"Visible",false,"FontName","Times New Roman","FontSize",20)
ax = gcf;
exportgraphics(ax,'./exported_images/colorbar_top.emf')
