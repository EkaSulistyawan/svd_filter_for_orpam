%% plot
titlelist = [['Raw'],"BPF","BPF+SVD wU", "BPF+SVD wV",...
    "Sym6-SURE","Sym4-DWT","Sym4-MODWT"];
figure("Position",[2012,113,1623,916])
tl = tiledlayout(4,7,'TileSpacing','tight');

seldat.raw = space3D;
seldat.bpf = space3Dbpf;
seldat.wu = denoisedSVDwU;
seldat.wv = denoisedSVDwV;
seldat.sym6 = denoisedPaper1Wavelet;
seldat.sym4dwt = denoisedPaper2DWT;
seldat.sym4modt = denoisedPaper2MODWT;

varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    hold on
    % view roi
    % xbox = [sgx(1), sgx(end), sgx(end), sgx(1)  , sgx(1)];
    % ybox = [sgy(1), sgy(1)  , sgy(end), sgy(end), sgy(1)];
    % plot(xbox,ybox,'Color','black','LineWidth',1.5); hold on
    % xbox = [bgx(1), bgx(end), bgx(end), bgx(1)  , bgx(1)];
    % ybox = [bgy(1), bgy(1)  , bgy(end), bgy(end), bgy(1)];
    % plot(xbox,ybox,'Color','green','LineWidth',1.5); hold off

    line([1 200],[9 9],'LineWidth',2,'Color','g','LineStyle','--')
    % line([68 165],[9 9],'LineWidth',2,'Color','c','LineStyle','--')
    line([10 30],[190 190],'LineWidth',2,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+14)
xaxis = (1:98)*250 / 1000;
plot(xaxis,normalize(im(68:165,9),'range'),'LineWidth',2)
axis tight
% ylim([0.05 0.8])
% xlabel(['Length (' char(181) 'm)'])
xlabel('Length ($\mu$m)','Interpreter','latex');
ylabel(' Norm. Intensity (a.u.)')
set(gca,'FontName','Times New Roman','FontSize',12)


nexttile(i+7)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:120)*250 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
% xlabel(['Length (' char(181) 'm)'])
% ylabel(['Depth (' char(181) 'm)'])
xlabel('Length ($\mu$m)','Interpreter','latex');
ylabel('Depth ($\mu$m)','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+21)
taxis = (1:1024)*1e6/Fs;
dat3d = seldat.(varnames{i});
sgAline_posX = 65;
sgAline_posY = 54;
bgAline_posX = 101;
bgAline_posY = 55;
plot(taxis,dat3d(:,sgAline_posX,sgAline_posY),'DisplayName','Signal','Color','Red','LineWidth',2);hold on;
plot(taxis,dat3d(:,bgAline_posX,bgAline_posY),'DisplayName','Unwanted','Color','Black','LineWidth',2);hold off
xlabel("Time ($\mu$s)",'Interpreter','latex')
ylabel("Intensity (a.u.)")
ylim([-6.5e-4 7.5e-4]);
set(gca,'FontName','Times New Roman','FontSize',12)


end

    p1 = [0.1 0.1];
    p2 = [0.15,0.8];
    annotation('textarrow',p1,p2,'String','Artifact','Color','cyan')

% cbh = colorbar;
% cbh.Layout.Tile = 'north';
% put legend
leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
leg.Layout.Tile = 'south';
ax = gcf;
exportgraphics(ax,'./exported_images/supp_rbc_artifact.emf')

%% custom
titlelist = [['Raw'],"BPF","BPF+SVD wU", "BPF+SVD wV",...
    "Sym6-SURE","Sym4-DWT","Sym4-MODWT"];
figure("Position",[2012,113,1623,916])
tl = tiledlayout(4,7,'TileSpacing','none');

seldat.raw = space3D;
seldat.bpf = space3Dbpf;
seldat.wu = denoisedSVDwU;
seldat.wv = denoisedSVDwV;
seldat.sym6 = denoisedPaper1Wavelet;
seldat.sym4dwt = denoisedPaper2DWT;
seldat.sym4modt = denoisedPaper2MODWT;

varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    hold on
    % view roi
    % xbox = [sgx(1), sgx(end), sgx(end), sgx(1)  , sgx(1)];
    % ybox = [sgy(1), sgy(1)  , sgy(end), sgy(end), sgy(1)];
    % plot(xbox,ybox,'Color','black','LineWidth',1.5); hold on
    % xbox = [bgx(1), bgx(end), bgx(end), bgx(1)  , bgx(1)];
    % ybox = [bgy(1), bgy(1)  , bgy(end), bgy(end), bgy(1)];
    % plot(xbox,ybox,'Color','green','LineWidth',1.5); hold off

    line([1 200],[9 9],'LineWidth',2,'Color','g','LineStyle','--')
    % line([68 165],[9 9],'LineWidth',2,'Color','c','LineStyle','--')
    line([10 30],[190 190],'LineWidth',2,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+14)
xaxis = (1:98)*250 / 1000;
plot(xaxis,im(68:165,9),'LineWidth',2)
axis tight
% ylim([0.05 0.8])
% xlabel(['Length (' char(181) 'm)'])
xlabel('Length ($\mu$m)','Interpreter','latex');
ylabel(' Norm. Intensity (a.u.)')
set(gca,'FontName','Times New Roman','FontSize',12)


nexttile(i+7)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:120)*250 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
if (i > 1)
    axis off
end
% xlabel(['Length (' char(181) 'm)'])
% ylabel(['Depth (' char(181) 'm)'])
xlabel('Length ($\mu$m)','Interpreter','latex');
ylabel('Depth ($\mu$m)','Interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',12)

nexttile(i+21)
taxis = (1:1024)*1e6/Fs;
dat3d = seldat.(varnames{i});
sgAline_posX = 65;
sgAline_posY = 54;
bgAline_posX = 101;
bgAline_posY = 55;
plot(taxis,dat3d(:,sgAline_posX,sgAline_posY),'DisplayName','Signal','Color','Red','LineWidth',2);hold on;
plot(taxis,dat3d(:,bgAline_posX,bgAline_posY),'DisplayName','Unwanted','Color','Black','LineWidth',2);hold off
xlabel("Time ($\mu$s)",'Interpreter','latex')
ylabel("Intensity (a.u.)")
ylim([-6.5e-4 7.5e-4]);
set(gca,'FontName','Times New Roman','FontSize',12)


end


% cbh = colorbar;
% cbh.Layout.Tile = 'north';
% put legend
leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
leg.Layout.Tile = 'south';
%% only two line
titlelist = [['Raw'],"BPF","BPF+SVD wU", "BPF+SVD wV",...
    "Sym6-SURE","Sym4-DWT","Sym4-MODWT"];
figure("Position",[2012,113,1623,916])
tl = tiledlayout(2,7,'TileSpacing','tight');

seldat.raw = space3D;
seldat.bpf = space3Dbpf;
seldat.wu = denoisedSVDwU;
seldat.wv = denoisedSVDwV;
seldat.sym6 = denoisedPaper1Wavelet;
seldat.sym4dwt = denoisedPaper2DWT;
seldat.sym4modt = denoisedPaper2MODWT;

varnames = fieldnames(seldat);


for i=1:numel(varnames)
nexttile(i)
im = cmodefunc(seldat.(varnames{i}));
imagesc(im');title(titlelist(i));colormap hot;axis off
if(i == 1)
    hold on
    % view roi
    xbox = [sgx(1), sgx(end), sgx(end), sgx(1)  , sgx(1)];
    ybox = [sgy(1), sgy(1)  , sgy(end), sgy(end), sgy(1)];
    plot(xbox,ybox,'Color','black','LineWidth',1.5); hold on
    xbox = [bgx(1), bgx(end), bgx(end), bgx(1)  , bgx(1)];
    ybox = [bgy(1), bgy(1)  , bgy(end), bgy(end), bgy(1)];
    plot(xbox,ybox,'Color','green','LineWidth',1.5); hold off
    line([1 200],[9 9],'LineWidth',2,'Color','g','LineStyle','--')
    % line([68 165],[9 9],'LineWidth',2,'Color','c','LineStyle','--')
    line([10 30],[190 190],'LineWidth',2,'Color','white')
end
set(gca,'FontName','Times New Roman','FontSize',12)


nexttile(i+7)
im = bmodefunc(seldat.(varnames{i}));
xaxis = (1:120)*250 / 1000; 
yaxis = (1:1024)*1e6*1480/Fs;
imagesc(xaxis,yaxis,im);colormap hot;
xlabel(['Length (' char(181) 'm)'])
ylabel(['Depth (' char(181) 'm)'])
set(gca,'FontName','Times New Roman','FontSize',12)



end


% cbh = colorbar;
% cbh.Layout.Tile = 'north';
% put legend
% leg = legend('Orientation', 'Horizontal','FontName','Times New Roman','FontSize',12);
% leg.Layout.Tile = 'south';
% ax = gcf;
% exportgraphics(ax,'./exported_images/supp_rbc_artifact.emf')


%% view roi
im = cmodefunc(space3D);
imagesc(im');hold on
xbox = [sgx(1), sgx(end), sgx(end), sgx(1)  , sgx(1)];
ybox = [sgy(1), sgy(1)  , sgy(end), sgy(end), sgy(1)];
plot(xbox,ybox,'Color','black','LineWidth',1.5); hold on
xbox = [bgx(1), bgx(end), bgx(end), bgx(1)  , bgx(1)];
ybox = [bgy(1), bgy(1)  , bgy(end), bgy(end), bgy(1)];
plot(xbox,ybox,'Color','green','LineWidth',1.5); hold on
colormap hot;axis square
%% show image collection

load("rbc_artifact_USV.mat",'u','v','s','ubpf','vbpf','sbpf','sz')
%%
up = u; vp = v; sp = s;

% recalculate weight
[wu,wv,ws] = getweight(up,vp,sp);

%no filter 
figure('Position',[59,65,1061,914]);
nr = 6;nc=6;
tiledlayout(nr,nc,'TileSpacing','tight','Padding','none');

for ii=1:nr*nc
    nexttile
    imagesc(reshape(vp(:,ii),sz,sz)')
    axis off
    colormap hot

    hold on

    sg = up(:,ii);
    sg = normalize(sg,'range')*sz;
    %sg = filtfilt(b,a,sg);
    sg = resample(sg,sz,1024);
    plot(sg,'Color','Cyan','LineWidth',2);
    
    hold on
    bgc = [1 1 1 0.75];
    text(5,30,sprintf("#%d",ii),'FontSize',18,'BackgroundColor',bgc);
    text(110,30,sprintf("$\\Sigma=%.2f$",ws(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,130,sprintf("$w_U=%.2f$",wu(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,180,sprintf("$w_V=%.2f$",wv(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    
    
    
    hold off
end

%% filter 
upf = ubpf; vpf = vbpf; spf = sbpf;
[wu,wv,ws] = getweight(upf,vpf,spf);

figure('Position',[59,65,1061,914]);
nr = 6;nc=6;
tiledlayout(nr,nc,'TileSpacing','tight','Padding','none');

for ii=1:nr*nc
    nexttile
    imagesc(reshape(vpf(:,ii),sz,sz)')
    axis off
    colormap hot

    hold on

    sg = upf(:,ii);
    sg = normalize(sg,'range')*sz;
    %sg = filtfilt(b,a,sg);
    sg = resample(sg,sz,1024);
    plot(sg,'Color','Cyan','LineWidth',2);
    
    hold on
    bgc = [1 1 1 0.75];
    text(5,30,sprintf("#%d",ii),'FontSize',18,'BackgroundColor',bgc);
    text(110,30,sprintf("$\\Sigma=%.2f$",ws(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,130,sprintf("$w_U=%.2f$",wu(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    text(5,180,sprintf("$w_V=%.2f$",wv(ii)),'FontSize',12,'Interpreter','latex','BackgroundColor',bgc);
    
    
    hold off
end
%%
function p = cmodefunc(x)
    p = squeeze(max(abs(hilbert(x))));
    % p = (p - min(p(:))) / range(p(:));
    % p(p > 0.45) = 0.45;
    % p = (p - min(p(:))) / range(p(:));
end

function p = bmodefunc(x)
    hb = abs(hilbert(x));
    p = squeeze(hb(:,:,9));
    p = (p - min(p(:))) / range(p(:));
end

function [wU,wV,wS] = getweight(u,v,s)
    imsize = round(sqrt(size(v,1)));

    uu = normalize(u,1,"range")-0.5;
    uu = sum(uu.^2);
    wU = normalize(-uu,'range');
    %wU = wU' .* diag(s);
    
    % wV
    v3d = reshape(v,imsize,imsize,1024);
    ss = pagesvd(v3d,'econ');
    vv = sum(squeeze(ss));
    wV = normalize(-vv,'range');
    %wV = wV' .* diag(s);

    wS = diag(s);
end