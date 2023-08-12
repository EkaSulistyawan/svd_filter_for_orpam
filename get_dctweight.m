function [y] = get_dctweight(d)
% d is a 3D data
datsize = size(d,3);
y = zeros(1,datsize);
for i=1:datsize
    roi = squeeze(d(:,:,i));
    dctroi = dct2(abs(roi));
    dctroi(1,1) = 0;
    maxdctroi = max(dctroi(:));
    y(i) = maxdctroi;
end

