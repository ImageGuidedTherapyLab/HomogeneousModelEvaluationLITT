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
mu_eff_tag = [0 900];

[opt, LOOCV, fig_labels, total] = master_LOOCV_onlySS( datafilename, DSC_thresholds, mu_thresholds,naive_mu(1), mu_eff_tag);

survival_plot_onlySS (LOOCV.dice.values, LOOCV.dice.naive.val, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_mu);