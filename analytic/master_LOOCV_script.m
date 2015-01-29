close all
clear

opttype = 'bestfit50';          % Name the opttype
datafilename = 'data_summary_GPU.txt';
%datafilename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file
DSC_thresholds = [0.3 ];
DSC_thresholds = sort(DSC_thresholds);
mu_thresholds = [];  % Intervals: 0 to 100, 100 to 150, 150 to 200, 200
mu_thresholds = sort(mu_thresholds);


[opt, LOOCV, fig_labels] = master_LOOCV_onlySS ( opttype, datafilename, DSC_thresholds, mu_thresholds);

survival_plot (LOOCV.dice.values, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2);