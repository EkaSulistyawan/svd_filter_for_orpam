clear all
load("../../SharingPoint/data/cellular/compiled/17122022_USAF.mat")
space3D = space3D(:,1:200,1:200);
Fs = 5e9;

organized_show_image(space3D)

%% 
sg = space3D(:,103,28);
ns = space3D(:,49,100);

%%

L = sqrt(size(sg,1));
A=double(sg); 
B=double(ns); 
A = (A - min(A)) / range(A);
B = (B - min(B)) / range(B);
%B = A;
bins = linspace(0,1,L+1);
na = zeros(1,L);
nb = zeros(1,L);
nab = zeros(L,L);

for i=1:numel(bins)-1
    condA = (A > bins(i) & (A <= bins(i+1)));
    condB = (B > bins(i) & (B <= bins(i+1)));
    na(i) = sum(condA);
    nb(i) = sum(condB);
    for j=1:numel(bins)-1
        condA = (A > bins(i) & (A <= bins(i+1)));
        condB = (B > bins(j) & (B <= bins(j+1)));
        nab(i,j) = sum(condA.*condB);
    end
end
na = na / sum(na);
nb = nb / sum(nb);
nab = nab / sum(nab,'all');
na(na==0) = eps;
nb(nb==0) = eps;
nab(nab==0) = eps;
miVal = sum(nab.*log10(nab ./ (na'*nb)),'all');


figure;
tiledlayout(1,2);
nexttile()
imagesc(na'*nb);axis image;colormap jet
nexttile()
imagesc(nab);axis image;colormap jet

a = colorbar;
a.Layout.Tile = 'south';