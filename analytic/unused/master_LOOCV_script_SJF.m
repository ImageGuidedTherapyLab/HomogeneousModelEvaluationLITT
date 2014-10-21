% My script for master_LOOCV

opttype = 'bestfit50';
datafilename = 'datasummaryL2_10sourceNewton50.txt';
DSC_thresholds = [0.3 0.5 0.7];
mu_thresholds = [100 150 200];
Matlab_flag = 1;

if Matlab_flag ~=1 && Matlab_flag ~=0
    disp('Invalid Matlab_flag. Only 0 or 1 allowed')
    break
end

[stats, passes] = master_LOOCV ( opttype, datafilename, DSC_thresholds, mu_thresholds, Matlab_flag);