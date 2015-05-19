clear ii jj
close all


for jj=1:size(total_all(:,5),1)
    ii = index(jj);
    figure(ii);
    [AX h1 h2] = plotyy(total_all{ii,2}(:,1),total_all{ii,2}(:,2),total_all{ii,2}(:,1),total_all{ii,3}(:,7));
    legend([h1;h2],'L_2','DSC');
end