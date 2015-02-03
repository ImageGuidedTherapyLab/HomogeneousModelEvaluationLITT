function [total,best,datasummary] = arrange_total_dataset ( data_filename, mu_eff_tag );
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

cd ../../../MATLAB/Tests/direct_search
load total

cd ( path22 );

% ix=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1);
% total(ix,:) = [];
% 
% ix=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1);
% total(ix,:) = [];
% 
% ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
% total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0409'))==1);
total(ix,:) = [];

total_all = total;
if mu_eff_tag(1) == 1    % Eliminate the higher values
    [~, mx] = min( abs( total{1,2}(:,1) - mu_eff_tag(2))); 
    for ii=1:size(total,1)
        aa_iter = total{ii,2};
        aa_iter( mx:end,: ) = []; 
        total{ii,2} = aa_iter;
        aa_iter = total{ii,3};
        aa_iter( mx:end, : ) = [];
        total{ii,3} = aa_iter;
        [total{ii,4}(1), total{ii,4}(3) ] = min(total{ii,2}(:,2));
        total{ii,4}(2) = total{ii,2}( total{ii,4}(3),1);
        [total{ii,5}(1), total{ii,5}(3) ] = max(total{ii,3}(:,7));
        total{ii,5}(2) = total{ii,2}( total{ii,5}(3),1);
    end
    
elseif mu_eff_tag(1) ==2   % Eliminate the lower values
    [~, mx] = min( abs( total{1,2}(:,1) - mu_eff_tag(2))); 
    for ii=1:size(total,1)
        aa_iter = total{ii,2};
        aa_iter( 1:mx, : ) = []; 
        total{ii,2} = aa_iter;
        aa_iter = total{ii,3};
        aa_iter( 1:mx, : ) = [];
        total{ii,3} = aa_iter;
        [total{ii,4}(1), total{ii,4}(3) ] = min(total{ii,2}(:,2));
        total{ii,4}(2) = total{ii,2}( total{ii,4}(3),1);
        [total{ii,5}(1), total{ii,5}(3) ] = max(total{ii,3}(:,7));
        total{ii,5}(2) = total{ii,2}( total{ii,5}(3),1);
    end
end

num_runs = size(total_all,1);
num_mu   = size(total_all{1,2},1);
L2 = zeros( num_mu, num_runs);
dice = L2;
for ii = 1:num_runs
    
    L2(:,ii) = total_all{ii,2}(:,2);
    dice(:,ii) = total_all{ii,3}(:,7);

end
clear ii

L2_mean = mean( L2,2);
L2_median = median( L2,2);
L2_stDev  = std ( L2')';
dice_mean = mean( dice,2);
dice_median = median( dice,2);
dice_stDev  = std ( dice')';

[best.all.L2_mean.L2 best.all.L2_mean.ix]= min( L2_mean);
[best.all.L2_median.L2 best.all.L2_median.ix] = min( L2_median );
[best.all.dice_mean.DSC best.all.dice_mean.ix] = max( dice_mean );
[best.all.dice_median.DSC best.all.dice_median.ix] = max( dice_median );

figure; h_title = title( 'L_2 mean, median, and st. dev.'); hold all;
[h1] = plot (total_all{1,2}(:,1), [L2_mean L2_median L2_stDev]);
legend( h1, 'L_2 mean', 'L_2 median', 'L_2 st. dev.'); hold off;

figure; h_title = title( 'DSC mean, median, and st. dev.'); hold all;
[h1] = plot (total_all{1,2}(:,1), [dice_mean dice_median dice_stDev]);
legend( h1, 'DSC mean', 'DSC median', 'DSC st. dev.'); hold off;

best.all.L2_mean.val = zeros( num_runs,1);
best.all.L2_median.val = best.all.L2_mean.val;
best.all.dice_mean.val = best.all.L2_mean.val;
best.all.dice_median.val = best.all.L2_mean.val;
for ii = 1:num_runs
    
    best.all.L2_mean.val(ii)        = total_all{ii,3}(best.all.L2_mean.ix,    7);
    best.all.L2_median.val(ii)      = total_all{ii,3}(best.all.L2_median.ix,  7);
    best.all.dice_mean.val(ii)      = total_all{ii,3}(best.all.dice_mean.ix,  7);
    best.all.dice_median.val(ii)    = total_all{ii,3}(best.all.dice_median.ix,7);
    
end


end  