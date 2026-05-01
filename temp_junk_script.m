load('D:\2026-01-16_MC_SC_17\post\masks\s03L\s03L_masks_QC.mat')

z_span = masks_table.("Z-span");
Result = masks_table.("Result");
uniq = unique(z_span);
frac_good = zeros(1,length(uniq));

for i=1:length(uniq)
    mask = (z_span == uniq(i));
    result_uniq = Result(mask);
    num_uniq = sum(mask);
    num_good = sum(strcmp(result_uniq, 'good'));   
    frac_good(i) = num_good/num_uniq;
end

x = uniq;
y = frac_good;

bar(x, y)
xlabel('Z-span (planes)')
ylabel('Fraction of masks passed')
set(gca, 'FontSize', 20)