close all
clear
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
choice = 4;  % 1 = mu; 2 = perf; 3 = cond;
toss_choice = 1;

if choice == 1;   % mu
    naive_var = [180 1];
    var_thresholds = [ ];  %  Eg: [ 100 150 200 ]; Intervals: 0 to 100, 100 to 150, 150 to 200, 200
    var_tag = [1 400];  % Splits the data at this proscribed point (index 2). Index 1 is 1 = save low side, 2 = save high side, 0 = keep all
    
elseif choice ==2;   % perf
    naive_var = [6 1];
    var_thresholds = [ ];  %  Eg: [ 100 150 200 ]; Intervals: 0 to 100, 100 to 150, 150 to 200, 200
    var_tag = [0 6];  % Splits the data at this proscribed point (index 2). Index 1 is 1 = save low side, 2 = save high side, 0 = keep all
    
    
elseif choice ==3;   % cond
    naive_var = [0.527 1];
    var_thresholds = [ ];  %  Eg: [ 100 150 200 ]; Intervals: 0 to 100, 100 to 150, 150 to 200, 200
    var_tag = [0 0.7];  % Splits the data at this proscribed point (index 2). Index 1 is 1 = save low side, 2 = save high side, 0 = keep all
    
elseif choice ==4;   % cond
    naive_var = [0.527 1];
    var_threholds = [];
    var_tag = [300 0; 0 0];
    %var_threholds = [300, 1; 0 0];
    
elseif choice ==5;   % mu + perf random
    naive_var = [0.527 1];
    var_threholds = [];
    var_tag = [300 0; 0 0];
    %var_threholds = [300, 1; 0 0];
  
end

%opttype = 'bestfit50';          % Name the opttype
datafilename = 'data_summary_GPU_L2.txt';
%datafilename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file
DSC_thresholds = [ -1 ];
DSC_thresholds = sort(DSC_thresholds);
opt_tag = 1; % DSC is 1; L2 is 2


if choice == 5||choice ==4
    [total,total_all,summary,toss_ix] = arrange_total_dataset_rand ( datafilename,var_tag, toss_choice );
    [opt, LOOCV, fig_labels] = master_LOOCV_onlySS_rand( total, DSC_thresholds,naive_var(1), opt_tag,choice);
    cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/libraries/troubleshoot_global_v_MC
    if toss_choice == 0
        aa = load ('mu_global_all');
    elseif toss_choice ==1
        aa= load('mu_global_21');
    end
    cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
    bb=hist2(opt.var.values{1}(:,1), opt.var.values{1}(:,2), 15);
    figure;imagesc(bb); xlabel('perf');ylabel('mu');
    survival_plot_onlySS_choice (LOOCV.dice.values, aa.LOOCV.dice.naive.val, LOOCV.var.run, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_var, choice);
    cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/libraries/troubleshoot_global_v_MC
    if toss_choice == 0
        
        if var_tag(1,2)==0
            save('MC_2para_all.mat', 'total','LOOCV','opt');
        else
            save ('MC_2para_all_no_var_thresh.mat','total','LOOCV','opt');
        end
        
    elseif toss_choice==1
        if var_tag(1,2)==0
            save('MC_2para_21_all.mat','total','LOOCV','opt');
        else
            save('MC_2para_21_no_var_thresh.mat','total','LOOCV','opt');
        end

    end
    cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
    
% elseif choice ==4
%     MC_figs_for_glory
    
    
    
else
    [total,total_all,summary] = arrange_total_dataset ( datafilename, var_tag, choice, toss_choice );
    [opt, LOOCV, fig_labels] = master_LOOCV_onlySS_choice( total, DSC_thresholds, var_thresholds,naive_var(1), opt_tag);
    survival_plot_onlySS_choice (LOOCV.dice.values, LOOCV.dice.naive.val, LOOCV.var.run, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_var, choice);
end


%[opt, LOOCV, fig_labels,best] = master_LOOCV_onlySS( total, DSC_thresholds, var_thresholds,naive_var(1), opt_tag, best, summary);



%[opt, LOOCV, fig_labels, total, total_all] = master_LOOCV_onlySS_naive_repeat( datafilename, DSC_thresholds, mu_thresholds,naive_mu(1), mu_eff_tag, opt_tag);

%naive_rerun( LOOCV.dice.values, LOOCV.dice.naive.val, LOOCV.mu_eff.run, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2, naive_mu