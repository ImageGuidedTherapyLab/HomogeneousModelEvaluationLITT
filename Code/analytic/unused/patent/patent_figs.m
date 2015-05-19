clear
close all
cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
load ('plot_pattern_amalgam_76_77.mat');

temp_paths = Study_paths25;
toss_index_LOOCV25 = find(dice_LOOCV25<0.7);
temp_paths(toss_index_LOOCV25,:) = [];
paths25_pass = temp_paths;
temp_dice = dice_LOOCV25;
temp_dice(toss_index_LOOCV25,:) = [];
dice25_pass = temp_dice;
temp_mueff = mu_eff25;
temp_mueff(toss_index_LOOCV25,:) = [];
mueff_pass = temp_mueff;

cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests/Patent_figs
handle1=figure(1); hist(mu_eff25);
print('-dpng','-r200','mu_eff_histSS.png')
handle2=figure(2); plot (thresholds, passes_LOOCV25,'LineWidth',5);
print('-dpng','-r200','dice_histSS.png')

num_good = length(dice25_pass);
pass_indices = zeros(num_good,1);
for ii = 1:num_good
    
    %pass_indices (ii) = find (ismember(Study_paths{:,2},Study_paths25{ii,2}));
    logics = strcmp(paths25_pass{ii,2},Study_paths25);
    pass_indices (ii) = find (logics(:,2));

end

[FOV_pass]=Patent_LOOCV_figures ( Study_paths25, mu_eff25, alpha25, best_iter25, opttype, Matlab_flag, pass_indices );