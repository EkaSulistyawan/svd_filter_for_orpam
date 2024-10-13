ref = load("\\10.32.224.8\Personal_folders\Eka\compiled\21032023_rbc_ref.mat").space3D;
ref = ref(:,1:200,1:200);
ref = awgn(ref,70);

Fs = 5e9;
[b,a] = butter(4,[10e6 100e6] / (Fs/2),'bandpass');
dat = filtfilt(b,a,ref);

%
[snrr,cnrr] = get_metrics(dat,'RBC');

%% get cmode
cmode_clean = squeeze(max(abs(hilbert(dat))));

%%
figure;
imagesc(cmode_clean)
axis image
colormap hot
colorbar
axis off

%%
sg_center_rbc = dat(:,127,3);
sg_peri_rbc = dat(:,113,12);
sg_bright = dat(:,137,112);
ns = dat(:,30,70);

% ftsg = abs(fft(sg_bright));
% ftns = abs(fft(ns));
% figure
% plot(10*log10(ftsg)); hold on
% plot(10*log10(ftns)); hold off
% xlim([0,100])
taxis = (1:1024) * 1e9 * (1/Fs);
plot(taxis,sg_bright);hold on
plot(taxis,ns);hold off
axis tight
%% implement svd
dat2d = reshape(dat,1024,200*200);
[u,s,v] = svd(dat2d,"econ");

%%
v3d = reshape(v',1024,200,200);

%%
% Create tiled layout with 2 rows and 3 columns
% tiledlayout(2, 3,"TileSpacing","none");

% Loop through idx from 1 to 6
for idx = 1:4
    % Select the appropriate tile for the current idx
    figure
    
    % Extract the 2D slice from 3D data
    vim = squeeze(v3d(idx,:,:));
    
    % Plot the image
    imagesc(vim)
    hold on
    
    % Get the signal
    sg = u(:,idx);
    %ori_range = max(sg) - min(sg);
    sg = resample(sg, 200, 1024);
    sg = normalize(-sg, 'range') * 200 ;
    
    % Overlay the signal on the image
    plot(sg, Color='cyan', LineWidth=5);
    
    % Set colormap, axis, and appearance
    colormap hot
    axis image off
end

%%
%% 
function status = view3Dpa(seldat,ths,ds)
    sz = size(seldat,2);
    dat2d = reshape(seldat,1024,sz*sz);
    dat2d = abs(hilbert(dat2d));
    dat1d = dat2d(:);
    
    dat1d = normalize(dat1d,'range');
    
    k = find(dat1d > ths);
    
    Fs = 5e9;
    kz = mod(k,1024)*1e9/Fs;
    kx = mod(ceil(k / 1024),sz)*200/1000;
    ky = ceil(ceil(k/1024)/sz)*200/1000;
    
    S = dat1d(k);
    C = dat1d(k);
    
    scatter3(kx(1:ds:end),ky(1:ds:end),kz(1:ds:end),S(1:ds:end).*10,C(1:ds:end));colormap hot
    colorbar

    status = true;
end

function [a,b] = get_metrics(D,datused)
    global_param;
    if datused == "Flat"
        sgx = sgx_flat2; 
        sgy = sgy_flat2;
        bgx = bgx_flat2; 
        bgy = bgy_flat2;
    elseif datused == "RBC"
        sgx = sgx_rbcIUS; 
        sgy = sgy_rbcIUS;
        bgx = bgx_rbcIUS; 
        bgy = bgy_rbcIUS;
    end
    
    sg = reshape(D(:,sgx,sgy),1024,size(sgx,2)*size(sgy,2));
    bg = reshape(D(:,bgx,bgy),1024,size(bgx,2)*size(bgy,2));
    
    taxis = (1:1024)*1e9/5e9;
    
    % plot mean & shade
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

    figure
    imagesc(sg)
    % plot mean & shade
    figure(Position=[438,602,765,214])
    alph = 0.3;
    plot_shade(taxis,sg(:,3),'Signal',[0.0 0.5 0.0],alph);hold on
    plot_shade(taxis,bg(:,1),'Noise',[1.0 0.0 0.5],alph);hold off
    axis tight
    set(gca,'FontName','Times','FontSize',25,'FontWeight','Bold')
    % ylim([-6 6]) % e-4
    ylim([-1 1]*10)
    xticklabels([])
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
        outp = outp + custom_snr2(a(:,i),b(:,i),i);
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
    maxdbns = max(ftb(idx));
    
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

function [outp] = custom_snr2(a,b,i)
plot(a)
    outp = snr(a(500:800),b(500:800));
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