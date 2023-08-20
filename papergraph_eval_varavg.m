load("var_avg_data/compiled_result.mat")

%% see the result for 1
cmodefunc = @(x)(squeeze(max(abs(hilbert(x)))));

cmoderaw = cmodefunc(result_avg_100.bpf);
maxval = max(cmoderaw,[],"all");
minval = min(cmoderaw,[],"all");
figure;imagesc(cmoderaw');axis image; colormap hot;caxis([minval maxval])
figure;imagesc(get_bmode(result_avg_100.bpf));colormap hot
%%
function p = get_bmode(x)
size(x)
    hb = abs(hilbert(x));
    p = squeeze(hb(:,60,:));
end