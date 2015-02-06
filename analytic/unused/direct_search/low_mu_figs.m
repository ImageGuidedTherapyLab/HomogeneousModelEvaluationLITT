clear ii jj
close all
aa = cell2mat( total_all(:,5));
index=find( (aa(:,3)<700)==1);

for jj=1:length(index)
    ii = index(jj);
    figure(ii);
    [AX h1 h2] = plotyy(total_all{ii,2}(:,1),total_all{ii,2}(:,2),total_all{ii,2}(:,1),total_all{ii,3}(:,7));
    legend([h1;h2],'L_2','DSC');
end