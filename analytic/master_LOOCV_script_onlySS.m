close all
clear
clc

opttype = 'bestfit50';          % Name the opttype
datafilename = 'data_summary_GPU_L2.txt';
%datafilename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file
DSC_thresholds = [ 0.3 ];
DSC_thresholds = sort(DSC_thresholds);
mu_thresholds = [ ];  % Intervals: 0 to 100, 100 to 150, 150 to 200, 200
mu_thresholds = sort(mu_thresholds);
naive_mu = [ 180 1];
mu_eff_tag = [2 400];
opt_tag = 2; % DSC is 1; L2 is 2

[total,total_all] = arrange_total_dataset ( datafilename, DSC_thresholds, mu_thresholds,naive_mu(1), mu_eff_tag, opt_tag );

[opt, LOOCV, fig_labels, total, total_all] = master_LOOCV_onlySS( datafilename, DSC_thresholds, mu_thresholds,naive_mu(1), mu_eff_tag, opt_tag);

[best] = find_naive_guess (total_all);

survival_plot_onlySS (LOOCV.dice.values, LOOCV.dice.naive.val, LOOCV.mu_eff.run, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_mu, best);

[opt, LOOCV, fig_labels, total, total_all] = master_LOOCV_onlySS_naive_repeat( datafilename, DSC_thresholds, mu_thresholds,naive_mu(1), mu_eff_tag, opt_tag);

%naive_rerun( LOOCV.dice.values, LOOCV.dice.naive.val, LOOCV.mu_eff.run, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_mu