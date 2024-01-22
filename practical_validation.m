clear all
load("../../Tohoku/Research Project/CS/Matlab/data/cellular/compiled/17122022_USAF.mat")
sz = size(space3D,2)-1;
space3D = space3D(:,1:sz,1:sz);
cmodefunc = @(x)(squeeze(max(abs(hilbert(space3D)))));
figure;
imagesc(cmodefunc(space3D));
axis image

Fs = 5e9;

space2D = reshape(space3D,1024,sz^2);
%space2D = normalize(space2D,1,"zscore");

% % BPF
fcl = 10e6;
fch = 100e6;
[b,a] = butter(4, [fcl/Fs*2 fch/Fs*2]);
space2D2 = filtfilt(b,a,space2D);


[UP,SP,VP] = svd(space2D,"econ");
[UP2,SP2,VP2] = svd(space2D2,"econ");

%% build H
h = UP(:,1);
T=1024;
H = zeros(T,T);
for i=1:T
    H(i,:) = circshift(h,i);
end
H = circshift(H,512,2);

[UH,SH,VH] = svd(H);


h = UP2(:,1);
T=1024;
H = zeros(T,T);
for i=1:T
    H(i,:) = circshift(h,i);
end
H = circshift(H,512,2);

[UH2,SH2,VH2] = svd(H);

%%
figure
tiledlayout(2,2)
nexttile
imagesc(abs(fft(UP)));
nexttile
imagesc(abs(fft(UH)));
nexttile
imagesc(abs(fft(UP2)));
nexttile
imagesc(abs(fft(UH2)));

%%
SH = sort(abs(fft(UP(:,1))),'descend');
SH2 = sort(abs(fft(filtfilt(b,a,UP(:,1)))),'descend');

figure
semilogy(diag(SP),'DisplayName','$\Sigma_P$');hold on
semilogy((SH),'DisplayName',"$\Sigma_H$");hold on
semilogy(diag(SP2),'DisplayName','Filtered $\Sigma_P$');hold on
semilogy((SH2),'DisplayName',"Filtered $\Sigma_H$");hold off
legend('Interpreter','latex');
axis square tight

%%
figure;
nr = 4;
nc = 4;
tiledlayout(nr,nc)

for j=1:nr*nc
    i=j+0;
    nexttile
    imagesc(reshape(VP2(:,i),sz,sz));axis image
    hold on

    sg = UP2(:,i);
    sg = normalize(sg,'range') - 0.5;
    sg = sg*sz+sz/2;
    sg = resample(sg,sz,1024);
    plot(sg,'Color','red')
    hold off
end