clear
close all
cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
load ('plot_pattern_patent_amalgam.mat');

temp_paths = Study_paths3Hi;
toss_index_LOOCV3 = find(diceHi_LOOCV3<0.7);
temp_paths(toss_index_LOOCV3,:) = [];
paths3_pass = temp_paths;

temp_dice = diceHi_LOOCV3;
temp_dice(toss_index_LOOCV3,:) = [];
dice3_pass = temp_dice;

temp_mueff = mu_eff3Hi;
temp_mueff(toss_index_LOOCV3,:) = [];
mueff3_pass = temp_mueff;


temp_paths = Study_paths7Hi;
toss_index_LOOCV7 = find(diceHi_LOOCV7<0.7);
temp_paths(toss_index_LOOCV7,:) = [];
paths7_pass = temp_paths;

temp_dice = diceHi_LOOCV7;
temp_dice(toss_index_LOOCV7,:) = [];
dice7_pass = temp_dice;

temp_mueff = mu_eff7Hi;
temp_mueff(toss_index_LOOCV7,:) = [];
mueff7_pass = temp_mueff;

cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests/Patent_figs
handle(1)=figure(1); hist(mu_eff3Hi);
print('-dpng','-r200','mu_eff3_histSS.png')
handle(2)=figure(2); plot (thresholds, passes_LOOCV3Hi,'LineWidth',5);
print('-dpng','-r200','dice3_histSS.png')

handle(3)=figure(3); hist(mu_eff7Hi);
print('-dpng','-r200','mu_eff7_histSS.png')
handle(2)=figure(4); plot (thresholds, passes_LOOCV7Hi,'LineWidth',5);
print('-dpng','-r200','dice7_histSS.png')

num_good3 = length(dice3_pass);
num_good7 = length(dice7_pass);
pass_indices3 = zeros(num_good3,1);
pass_indices7 = zeros(num_good7,1);
for ii = 1:num_good7
    
    %pass_indices (ii) = find (ismember(Study_paths{:,2},Study_paths25{ii,2}));
    logics7 = strcmp(paths7_pass{ii,2},Study_paths7Hi);
    pass_indices7 (ii) = find (logics7(:,2));
    
    if ii <= num_good3
        logics3 = strcmp(paths3_pass{ii,2},Study_paths3Hi);
        pass_indices3 (ii) = find (logics3(:,2));

    end

end

[FOV_pass3]=Patent_LOOCV_figures ( Study_paths3Hi, mu_eff3Hi, alpha3Hi, best_iter3Hi, opttype2, Matlab_flag, pass_indices3, 0 );
[FOV_pass7]=Patent_LOOCV_figures ( Study_paths7Hi, mu_eff7Hi, alpha7Hi, best_iter7Hi, opttype2, Matlab_flag, pass_indices7, 1 );
