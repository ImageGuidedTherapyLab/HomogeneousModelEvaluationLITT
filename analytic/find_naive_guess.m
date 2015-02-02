function [best] = find_naive_guess ( total );

num_runs = size(total,1);
num_mu   = size(total{1,2},1);
L2 = zeros( num_mu, num_runs);
dice = L2;
for ii = 1:num_runs
    
    L2(:,ii) = total{ii,2}(:,2);
    dice(:,ii) = total{ii,3}(:,7);

end
clear ii

L2_mean = mean( L2,2);
L2_median = median( L2,2);
L2_stDev  = std ( L2')';
dice_mean = mean( dice,2);
dice_median = median( dice,2);
dice_stDev  = std ( dice')';

[best.L2_mean.val best.L2_mean.ix]= max( L2_mean);
[best.L2_median.val best.L2_median.ix] = max( L2_median );
[best.dice_mean.val best.dice_mean.ix] = max( dice_mean );
[best.dice_median.val best.dice_median.ix] = max( dice_median );

figure; h_title = title( 'L_2 mean, median, and st. dev.'); hold all;
[h1] = plot (total{1,2}(:,1), [L2_mean L2_median L2_stDev]);
legend( h1, 'L_2 mean', 'L_2 median', 'L_2 st. dev.'); hold off;

figure; h_title = title( 'DSC mean, median, and st. dev.'); hold all;
[h1] = plot (total{1,2}(:,1), [dice_mean dice_median dice_stDev]);
legend( h1, 'DSC mean', 'DSC median', 'DSC st. dev.'); hold off;
% [ax, h1, h2] = plotyy ( total{1,2}(:,1), [L2_mean L2_median], total{1,2}(:,1), L2_stDev, 'plot','plot' );
% legend( [h1;h2], 'L_2 mean', 'L_2 median', 'L_2 st. dev.'); hold off;

% figure; h_title = title( 'DSC mean, median, and st. dev.'); hold all;
% [ax, h1, h2] = plotyy ( total{1,2}(:,1), [dice_mean dice_median], total{1,2}(:,1), dice_stDev, 'plot','plot' );
% legend( [h1;h2], 'DSC mean', 'DSC median', 'DSC st. dev.'); hold off;

% plot ( total{1,2}(:,1), L2_mean  );
% legend_string = ['L2_mean'];
% legend( legend_string, 'Location','southwest');
% hold all;
% legend('-DynamicLegend', 'Location','southeast');
% legend_string = ['
% 
% plotyy ( total{1,2}(:,1), L2_median, total{1,2}(:,1), L2_stDev);
% hold off;
% 
% legend_string = strcat( ['Optimization']);
% legend( legend_string, 'Location','southwest');
% hold all
% legend('-DynamicLegend', 'Location','southwest');
% legend_string = strcat(['Naive '], num2str(naive_tag(1)), ' m^{-1}');
% plot (thresholds, pass_naive_opt, 'DisplayName',legend_string, 'LineWidth',5);
% hold off
% 
% figure; h_title = title( 'DSC mean, median, and st. dev.'); hold on;
% plot ( total{1,2}(:,1), dice_mean  );
% plotyy ( total{1,2}(:,1), dice_median, total{1,2}(:,1), dice_stDev);
% hold off;
