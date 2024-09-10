dat = load("\\10.32.224.8\Personal_folders\Eka\compiled\29112022_10umsphere.mat").space3D;
dat = dat(:,1:size(dat,2)-1,1:size(dat,2)-1);
imsz = size(dat,2);

Fs = 5e9;

fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
tic
space3Dbpf = filtfilt(b,a,dat);
tbpf = toc;

%%
cmode = squeeze(max(abs(hilbert(space3Dbpf))));
figure;imagesc((cmode));axis image off;colormap hot;



%% calculate depth map
envlp = abs(hilbert(space3Dbpf));
[maxval,Ibpf] = max(envlp,[],1);

%%
figure
plot((1:1024)*1e9/Fs,envlp(:,32,81))
%%

dmap = squeeze(Ibpf)'*1e9/Fs;
dmin = 70;
dmax = 85;
mask = (dmap < dmin) | (dmap > dmax);
dmap(mask) = 0;

jetmap = jet(256);
customjet = [0 0 0;jetmap];

figure;imagesc(dmap');
colormap(customjet);
axis image off
clim([dmin dmax])
cb = colorbar;
ylabel(cb,['Time Delay (ns)'])
%set(cb,'Ticks',[105, 115,125])
set(gca,'FontName','Times','FontSize',20,'FontWeight','bold')

%% plot the svd
dat2d = reshape(space3Dbpf,1024,imsz*imsz);
[u,s,v] = svd(dat2d,'econ');
v3d = reshape(v',1024,imsz,imsz);

%% plot the svd result
close all
idx = 6;

vim = squeeze(v3d(idx,:,:));

figure;
imagesc(vim)
hold on

% get the signal
sg = u(:,idx);
sg = resample(sg,imsz,1024);
sg = normalize(sg,'range') * imsz;
plot(sg,Color='cyan',LineWidth=5)
colormap hot
axis image off

sim = svd(vim,"econ","vector");
sim = sim';

%%
[snrr0ori,cnrr0ori] = get_metrics(space3Dbpf,'FLAT');
%% calculate SNR 

function [a,b] = get_metrics(D,datused)
    global_param;
    if datused == "FLAT"
        sgx = sgx_flat2; 
        sgy = sgy_flat2;
        bgx = bgx_flat2; 
        bgy = bgy_flat2;
    elseif datused == "RBC"
        sgx = sgx_rbcIUS; 
        sgy = sgy_rbcIUS;
        bgx = bgx_rbcIUS; 
        bgy = bgy_rbcIUS;
    elseif datused == "SPHERE"
        sgx = sgx_round2; 
        sgy = sgy_round2;
        bgx = bgx_round2; 
        bgy = bgy_round2;
    end
    
    sg = reshape(D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
    bg = reshape(D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
    
    taxis = (1:1024)*1e9/5e9;

    figure(Position=[438,602,765,214])
    alph = 0.3;
    plot_shade(taxis,squeeze(sg),'Signal',[0.0 0.5 0.0],alph);hold on
    plot_shade(taxis,squeeze(bg),'Noise',[1.0 0.0 0.5],alph);hold off
    axis tight
    set(gca,'FontName','Times','FontSize',25,'FontWeight','Bold')
    % ylim([-6 6]) % e-4
    ylim([-1 1]*10)
    %xticklabels([])
    yticklabels([])


    a = snr2d(sg,bg);
    b = get_cnr(D,sgx,sgy,bgx,bgy);
end

function [outp] = plot_shade(x,dat,nm,c1,ca)
    avg = mean(squeeze(dat),2)*1e4;
    stdv = std(squeeze(dat),[],2)*1e4;

    curve1 = avg + stdv;
    curve2 = avg - stdv;

    x2 = [x, fliplr(x)];
    inBetween = [curve1; flipud(curve2)];
    plot(x, avg,'Color', c1, 'LineWidth', 4,'DisplayName',strcat(['Avg.' nm]));hold on
    fill(x2, inBetween,c1,'FaceAlpha',ca,'EdgeColor','none','DisplayName',strcat(['Std.' nm]));hold on
    
end

function [outp] = snr2d(a,b)
    sz = size(a,2);
    outp = 0;
    for i =1:sz
        outp = outp + custom_snr(a(:,i),b(:,i),i);
    end
    outp = outp / sz;
end

function [a] = range(p)
    a = max(p,[],"all") - min(p,[],'all');
end


function [outp] = custom_snr(a,b,i)
    lowlim = 10e6;
    hghlim = 100e6;

    faxis = (1:512)*5e9/1024;

    idx = find((faxis > lowlim) & (faxis < hghlim));
    fta = abs(fft(a));
    ftb = abs(fft(b)); 
    maxdbsg = max(fta(idx));
    maxdbns = std(ftb(idx));
    
    % if i==1
    %     disp(i)
    %     figure;
    %     plot(fta(idx));hold on
    %     plot(ftb(idx));hold off
    %     title(i)
    % end
    outp = mag2db(maxdbsg/maxdbns);
    % fprintf("%.2f\n",outp);
end

function [outp] = snrv2(a,b,i)
    lowlim = 10e6;
    hghlim = 100e6;

    faxis = (1:512)*5e9/1024;

    idx = find((faxis > lowlim) & (faxis < hghlim));
    fta = abs(fft(a));
    ftb = abs(fft(b)); 
    maxdbsg = max(fta(idx));
    maxdbns = std(ftb(idx));
    
    outp = mag2db(maxdbsg/maxdbns);

end

function [a] = get_cnr(D,sgx,sgy,bgx,bgy)

    hb = abs(hilbert(D));
    cmode = squeeze(max(hb));
    sgmu = mean(cmode(sgx,sgy),"all");
    bgmu = mean(cmode(bgx,bgy),"all");
    sgsd = std(cmode(sgx,sgy),0,"all");
    bgsd = std(cmode(bgx,bgy),0,"all");
    a = (sgmu - bgmu) / (sqrt((sgsd^2 + bgsd^2) / 2));
    a = mag2db(a);
end