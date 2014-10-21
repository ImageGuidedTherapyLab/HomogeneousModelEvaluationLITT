% My script for master_LOOCV
close all
clear

opttype = 'bestfit50';
datafilename = 'datasummaryL2_10sourceNewton50.txt';
DSC_thresholds = [0.3 0.5 0.7];
mu_thresholds = [100 150 200];
Matlab_flag = 1;

if Matlab_flag ~=1 && Matlab_flag ~=0
    disp('Invalid Matlab_flag. Only 0 or 1 allowed')
    break
end

[opt, LOOCV] = master_LOOCV ( opttype, datafilename, DSC_thresholds, mu_thresholds, Matlab_flag);

survival_plot (LOOCV.dice.values, opt.dice.all);

% thresholds = linspace ( 0.0, 1, 10001);
% %passes_LOOCV8 = zeros (10001,1);
% passes_LOOCV7 = zeros (10001,1);
% % passes_LOOCV25 = zeros (10001,1);
% passes_LOOCV3 = zeros(10001,1);
% % passes_LOOCV_Low7 = zeros(10001,1);
% % passes_LOOCV_High7 = zeros(10001,1);
% passes_all = zeros(10001,1);
% %passes_iter = passes_LOOCV8;
% for ii = 1:10001
%     %passes_iter   (ii) = sum ( dice_values_iter > thresholds(ii));
%     passes_LOOCV7 (ii) = sum ( dice_LOOCV7 > thresholds(ii));
% %     passes_LOOCV25 (ii) = sum ( dice_LOOCV25 > thresholds(ii));
%     passes_LOOCV3 (ii) = sum ( dice_LOOCV3 > thresholds(ii));
% %     passes_LOOCV_Low7 (ii) = sum ( dice_LOOCV_Low7 > thresholds(ii));
% %     passes_LOOCV_High7 (ii) = sum ( dice_LOOCV_High7 > thresholds(ii));
%     passes_all (ii) = sum ( datasummary(:,7) > thresholds(ii));
%     %passes_LOOCV8 (ii) = sum ( dice_LOOCV8 > thresholds(ii));
% end
% 
% %passes_iter_AUC = sum (passes_iter) ./ (10001 * stats_mu_iter.n) ;  % The AUC is actually the same as the mean
% passes_LOOCV7_AUC = sum (passes_LOOCV7) ./ (10001 * stats_mu7.n);
% % passes_LOOCV7_AUC_Low = sum (passes_LOOCV_Low7) ./ (10001 * stats_muLow7.n);
% % passes_LOOCV7_AUC_High = sum (passes_LOOCV_High7) ./ (10001 * stats_muHigh7.n);
% %passes_LOOCV8_AUC = sum (passes_LOOCV8) ./ (10001 * stats_mu8.n);
% %figure(1); plot (thresholds, passes_LOOCV8,'LineWidth',5);
% figure; plot (thresholds, passes_LOOCV7);
% % figure; plot (thresholds, passes_LOOCV25,'LineWidth',5);
% figure; plot (thresholds, passes_LOOCV3,'LineWidth',5);
% % figure; plot (thresholds, passes_LOOCV_Low7);
% % figure; plot (thresholds, passes_LOOCV_High7);
% figure; plot (thresholds, passes_all);
