load("saved_others\melanoma_test.mat");

%% quick function
bmodefunc = @(d,x)(abs(hilbert(squeeze(d(:,:,x)))));
cmodefunc = @(x)(squeeze(max(abs(hilbert(x)))));
alinefunc = @(d,x,y)(d(:,x,y));

bmode_pos = 80;
sgAline_posX = 40;
sgAline_posY = 74;
bgAline_posX = 94;
bgAline_posY = 25;

%% bmode slice
figure("Position",[8,63,1255,928])
tiledlayout(3,4,"TileSpacing","tight");

% cmode
nexttile
imagesc(cmodefunc(space3Dbpf)');axis off;colormap hot
line([1 125],[80 80],'LineWidth',2,'Color','g','LineStyle','--')
nexttile;
imagesc(cmodefunc(denoisedSVDgamma0)');axis off;colormap hot
annotation('arrow',[0.33 0.36],[0.73 0.76],'Color','c')
nexttile
imagesc(cmodefunc(denoisedSVDgammahalf)');axis off;colormap hot
nexttile
imagesc(cmodefunc(denoisedSVDgamma1)');axis off;colormap hot
line([100 120],[10 10],'LineWidth',2,'Color','white')
% colorbar = 20*400 nm

% bmode
xaxis = (1:125)*400 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;

nexttile
imagesc(xaxis,yaxis,bmodefunc(space3Dbpf,bmode_pos));colormap hot; 
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)
nexttile
imagesc(xaxis,yaxis,bmodefunc(denoisedSVDgamma0,bmode_pos));colormap hot; 
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)
nexttile
imagesc(xaxis,yaxis,bmodefunc(denoisedSVDgammahalf,bmode_pos));colormap hot; 
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)
nexttile
imagesc(xaxis,yaxis,bmodefunc(denoisedSVDgamma1,bmode_pos));colormap hot; 
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)

% Aline
taxis = (1:1024)*1e9/Fs;

nexttile
plot(taxis,alinefunc(space3Dbpf,sgAline_posX,sgAline_posY));hold on
plot(taxis,alinefunc(space3Dbpf,bgAline_posX,bgAline_posY));hold off
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([-2e-3 4e-3])
set(gca,'FontName','Times New Roman','FontSize',12)
nexttile
plot(taxis,alinefunc(denoisedSVDgamma0,sgAline_posX,sgAline_posY));hold on
plot(taxis,alinefunc(denoisedSVDgamma0,bgAline_posX,bgAline_posY));hold off
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([-2e-3 4e-3])
set(gca,'FontName','Times New Roman','FontSize',12)
nexttile
plot(taxis,alinefunc(denoisedSVDgammahalf,sgAline_posX,sgAline_posY));hold on
plot(taxis,alinefunc(denoisedSVDgammahalf,bgAline_posX,bgAline_posY));hold off
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([-2e-3 4e-3])
set(gca,'FontName','Times New Roman','FontSize',12)
nexttile
plot(taxis,alinefunc(denoisedSVDgamma1,sgAline_posX,sgAline_posY));hold on
plot(taxis,alinefunc(denoisedSVDgamma1,bgAline_posX,bgAline_posY));hold off
xlabel("Time (ns)")
ylabel("Intensity (a.u.)")
ylim([-2e-3 4e-3])
set(gca,'FontName','Times New Roman','FontSize',12)

axis tight

ax = gcf;
exportgraphics(ax,'./exported_images/fig_melanoma.emf')