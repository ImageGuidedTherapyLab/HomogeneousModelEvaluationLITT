close all
clear
clc

choice =1;  % 1 = mu; 2 = perf; 3 = cond;

if choice == 1;   % mu
    naive_var = [180 1];
    var_thresholds = [ ];  %  Eg: [ 100 150 200 ]; Intervals: 0 to 100, 100 to 150, 150 to 200, 200
    var_tag = [1 400];  % Splits the data at this proscribed point (index 2). Index 1 is 1 = save low side, 2 = save high side
    
elseif choice ==2;   % perf
    naive_var = [6 1];
    var_thresholds = [ ];  %  Eg: [ 100 150 200 ]; Intervals: 0 to 100, 100 to 150, 150 to 200, 200
    var_tag = [1 6];  % Splits the data at this proscribed point (index 2). Index 1 is 1 = save low side, 2 = save high side
    
    
elseif choice ==3;   % cond
    naive_var = [0.527 1];
    var_thresholds = [ ];  %  Eg: [ 100 150 200 ]; Intervals: 0 to 100, 100 to 150, 150 to 200, 200
    var_tag = [0 0.7];  % Splits the data at this proscribed point (index 2). Index 1 is 1 = save low side, 2 = save high side
    
end

var_thresholds = sort(var_thresholds);
opttype = 'bestfit50';          % Name the opttype
datafilename = 'data_summary_GPU_L2.txt';
%datafilename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file
DSC_thresholds = [ 0.3 ];
DSC_thresholds = sort(DSC_thresholds);
opt_tag = 1; % DSC is 1; L2 is 2

[total,total_all,summary] = arrange_total_dataset ( datafilename, var_tag, choice );
[opt, LOOCV, fig_labels] = master_LOOCV_onlySS_choice( total, DSC_thresholds, var_thresholds,naive_var(1), opt_tag, choice);



%[opt, LOOCV, fig_labels,best] = master_LOOCV_onlySS( total, DSC_thresholds, var_thresholds,naive_var(1), opt_tag, best, summary);

survival_plot_onlySS_choice (LOOCV.dice.values, LOOCV.dice.naive.val, LOOCV.var.run, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_var, choice);

%[opt, LOOCV, fig_labels, total, total_all] = master_LOOCV_onlySS_naive_repeat( datafilename, DSC_thresholds, mu_thresholds,naive_mu(1), mu_eff_tag, opt_tag);

%naive_rerun( LOOCV.dice.values, LOOCV.dice.naive.val, LOOCV.mu_eff.run, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_mu