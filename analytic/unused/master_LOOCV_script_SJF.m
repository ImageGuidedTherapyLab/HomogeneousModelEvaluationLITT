% My script for master_LOOCV
close all
clear

opttype = 'bestfit50';
datafilename = 'datasummaryL2_10sourceNewton50.txt';
% DSC_thresholds = [0.3 0.5 0.7];
% mu_thresholds = [100 150 200];
DSC_thresholds = [0.7];   % At least one DSC threshold is required
mu_thresholds = [175];   % Zero or more mu_threholds is possible
Matlab_flag = 1;

if Matlab_flag ~=1 && Matlab_flag ~=0
    disp('Invalid Matlab_flag. Only 0 or 1 allowed')
    break
end

[opt, LOOCV] = master_LOOCV ( opttype, datafilename, DSC_thresholds, mu_thresholds, Matlab_flag);

survival_plot (LOOCV.dice.values, opt.dice.all);