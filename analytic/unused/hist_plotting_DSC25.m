% This script is meant to see how many datasets have solutions that are
% good. The following questions need to be answered:
% 1. Does a good solution exist? This will be done by doing max(dice) >0.7
% 2. How far (in mu_eff value) is the max Dice from the min obj fxn?
%              I.e. (mu_eff@maxDice - mu_eff@minObjFxn)
% 3. Where does the optimizer go? Ie (mu_eff@Opt - mu_eff@minObjFxn) and
%           (mu_eff@Opt - mu_eff@maxDice)
% 4. How does the maxDice depend on isotherm choice (temperature).

% run run PlanningValidation directory

% path(path,'/workarea/fuentes/github/DakotaApplications/PlanningValidation/analytic/unused')
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
% load all_new_Tmp_files.mat
% 
% data.labels = total_path(:,1);
% 
% 
%
% num_runs = size(total_path,1);
% val_dice = zeros (num_runs,1);
% index_dice = val_dice;
% data.data = zeros (5001,6,num_runs);
% for ii=1:num_runs
%     data.data(:,:,ii) = total_path{ii,2};
% end

% Load the brute force sensitivity study
clear
close all

% I have the list of good runs, now I need to code it. 
% 
% The input: vector of paths to whatever interested datasets and vector of optimized mu_eff values.
% (I assume that the inverse problem is run before the LOOCV begins).
% 
% The process:
% Its iterative. Heres one iteration of leave-one-out cross validation (LOOCV). Average the optimized mu_eff for all but one ablation. 
% Use the averaged optimized mu_eff for the reserved ablation. 
% That means a dakota.in file will be written with the default constants and the averaged optimized mu_eff value.
% That dakota.in file wile be read by the BHTE code, the BHTE code runs, and a prediction is made. 
% Then, the prediction 57 deg_C isotherm is compared to the MRTI 57 dec_C isotherm using a Dice similarity coefficient (DSC). That is one iteration.
% Overall process iterates the LOOCV algorithm such that there are n possible iterations for n ablation datasets resulting in n DSC values. 
% The DSC values mean and variance are calculated. The mean and variance are sent to the t-test to find the hypotheses results.
% 
% The output: The output is a binary acceptance/rejection of the null and alternative hypotheses.

% This script finds the best mu_eff for the different studies.
opttype = 'bestfit' ;
Matlab_flag = 0; % 0 means use FEM kernel; 1 means use MATLAB for kernel




% read  best_fit optimization data and store mu_eff and alpha
%datasummary = dlmread('ex_datasummary.txt',',',1,0);
datasummary = dlmread('datasummary.txt',',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);

Study_paths = cell (num_studies,2);
for ii = 1:num_studies
    
    Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end

clear ii
indexC = strfind(Study_paths,'Study0035');
toss_index_phantom = find(not(cellfun('isempty',indexC)));
Study_paths(toss_index_phantom,:)=[];
datasummary(toss_index_phantom,:)=[];

if Matlab_flag == 0
    dice_raw = datasummary(:,7);

elseif Matlab_flag ==1
    dice_raw = 1 - datasummary(:,7);
    
else
    disp('Invalid Matlab_flag. Only 0 or 1 allowed')
    break
end
toss_index7 = find(dice_raw<0.7);
toss_index25= find(dice_raw<0.25);
% toss_index8 = find(dice_raw<0.8);
% toss_indexLow = find(datasummary(:,4)>300);
% toss_indexHigh = find(datasummary(:,4)<300);

num_studies = size(datasummary,1);
linenum=7;
quality=zeros(num_studies,1);
for ii = 1:num_studies
    
    fin = fopen(strcat('./workdir/',Study_paths{ii,1},'/',Study_paths{ii,2},'/opt/', 'setup.ini'));
    qual_data = textscan(fin, '%s', 1, 'delimiter', '\n', 'headerlines', linenum-1);
    qual_data1 = char(qual_data{1});
    quality(ii) = qual_data1(11);
    fclose(fin);
end
quality=char(quality);
clear ii

mu_eff7 = datasummary(:,4);
mu_eff7 (toss_index7) = [];
mu_eff25 = datasummary(:,4);
mu_eff25 (toss_index25) = [];
% mu_eff8 = datasummary(:,4);
% mu_eff8 (toss_index8) = [];
% mu_effLow = datasummary(:,4);
% mu_effLow(toss_indexLow) = [];

% toss_indexLow7 = find(mu_eff7>300);
% toss_indexHigh7 = find(mu_eff7<300);
% mu_effLow7 = mu_eff7;
% mu_effLow7(toss_indexLow7) = [];
% mu_effHigh7 = mu_eff7;
% mu_effHigh7(toss_indexHigh7)=[];


alpha7 = datasummary(:,5);
alpha7 (toss_index7) = [];
alpha25 = datasummary(:,5);
alpha25 (toss_index25) = [];
% alpha8 = datasummary(:,5);
% alpha8 (toss_index8) = [];
% alphaLow7 = alpha7;
% alphaLow7(toss_indexLow7)=[];
% alphaHigh7=alpha7;
% alphaHigh7(toss_indexHigh7)=[];

best_iter7 = datasummary(:,3);
best_iter7 (toss_index7) = [];
best_iter25 = datasummary(:,3);
best_iter25 (toss_index25) = [];
% best_iter8 = datasummary(:,3);
% best_iter8 (toss_index8) = [];
% best_iterLow7 = best_iter7;
% best_iterLow7(toss_indexLow7)=[];
% best_iterHigh7=best_iter7;
% best_iterHigh7(toss_indexHigh7)=[];

dice7 = dice_raw;
dice7(toss_index7) = [];
dice25 = dice_raw;
dice25(toss_index25) = [];
% dice8 = dice_raw;
% dice8(toss_index8) = [];
% diceLow7 = dice7;
% diceLow7(toss_indexLow7)=[];
% diceHigh7=dice7;
% diceHigh7(toss_indexHigh7)=[];

quality7 = quality;
quality7(toss_index7) = [];
quality25 = quality;
quality25(toss_index25) = [];

stats_mu_raw = Descriptive_statistics( datasummary(:,4) );
stats_alpha_raw = Descriptive_statistics( datasummary(:, 5) );
stats_mu7 = Descriptive_statistics(mu_eff7);
stats_mu25 = Descriptive_statistics(mu_eff25);
% stats_muLow7 = Descriptive_statistics(mu_effLow7);
% stats_muHigh7 = Descriptive_statistics(mu_effHigh7);
stats_alpha7 = Descriptive_statistics(alpha7);
stats_alpha25 = Descriptive_statistics(alpha25);
% stats_mu8 = Descriptive_statistics(mu_eff8);
% stats_alpha8 = Descriptive_statistics(alpha8);
stats_dice_raw = Descriptive_statistics(dice_raw);
stats_dice7 = Descriptive_statistics(dice7);
stats_dice25 = Descriptive_statistics(dice25);
% stats_diceLow7 = Descriptive_statistics(diceLow7);
% stats_diceHigh7 = Descriptive_statistics(diceHigh7);


figure; hist(datasummary(:,4));
figure; hist(mu_eff7);
figure; hist(mu_eff25);
% figure; hist(mu_effLow7);
% figure; hist(mu_effHigh7);

% Remove Study_paths indices
temp_paths = Study_paths;
temp_paths(toss_index7,:) = [];
Study_paths7 = temp_paths;
temp_paths = Study_paths;
temp_paths(toss_index25,:) = [];
Study_paths25 = temp_paths;
% temp_paths = Study_paths7;
% temp_paths(toss_indexLow8,:)=[];
% Study_pathsLow8=temp_paths;
% temp_paths = Study_paths7;
% temp_paths(toss_indexLow7,:)=[];
% Study_pathsLow7=temp_paths;
% temp_paths = Study_paths7;
% temp_paths(toss_indexHigh7,:)=[];
% Study_pathsHigh7=temp_paths;
% temp_paths = Study_paths;
% temp_paths(toss_index8,:) = [];
% Study_paths8 = temp_paths;

%[ hh_raw, dice_values_LOOCV ] = LOOCV_t_test_DiceTemp ( Study_paths, datasummary(:,4) ,datasummary(:,5) , datasummary(:,3), opttype, Matlab_flag );
%[ hh7, dice_LOOCV7 ] = LOOCV_t_test_DiceTemp ( Study_paths7, mu_eff7, alpha7, best_iter7, opttype, Matlab_flag );
[ hh25, dice_LOOCV25 ] = LOOCV_t_test_DiceTemp ( Study_paths25, mu_eff25, alpha25, best_iter25, opttype, Matlab_flag );
% [ hhLow7, dice_LOOCV_Low7 ] = LOOCV_t_test_DiceTemp ( Study_pathsLow7, mu_effLow7, alphaLow7, best_iterLow7, opttype, Matlab_flag );
% [ hhHigh7, dice_LOOCV_High7 ] = LOOCV_t_test_DiceTemp ( Study_pathsHigh7, mu_effHigh7, alphaHigh7, best_iterHigh7, opttype, Matlab_flag );
%[ hh8, dice_LOOCV8 ] = LOOCV_t_test_DiceTemp ( Study_paths8, mu_eff8, alpha8, best_iter8, opttype, Matlab_flag );
%dice_LOOCV8_stats = Descriptive_statistics(dice_LOOCV8);
%dice_LOOCV7_stats = Descriptive_statistics(dice_LOOCV7);
% dice_LOOCV_Low7stats = Descriptive_statistics(dice_LOOCV_Low7);
% dice_LOOCV_High7stats = Descriptive_statistics(dice_LOOCV_High7);


thresholds = linspace ( 0.0, 1, 10001);
%passes_LOOCV8 = zeros (10001,1);
%passes_LOOCV7 = zeros (10001,1);
passes_LOOCV25 = zeros (10001,1);
% passes_LOOCV_Low7 = zeros(10001,1);
% passes_LOOCV_High7 = zeros(10001,1);
passes_all = zeros(10001,1);
%passes_iter = passes_LOOCV8;
for ii = 1:10001
    %passes_iter   (ii) = sum ( dice_values_iter > thresholds(ii));
    %passes_LOOCV7 (ii) = sum ( dice_LOOCV7 > thresholds(ii));
    passes_LOOCV25 (ii) = sum ( dice_LOOCV25 > thresholds(ii));
%     passes_LOOCV_Low7 (ii) = sum ( dice_LOOCV_Low7 > thresholds(ii));
%     passes_LOOCV_High7 (ii) = sum ( dice_LOOCV_High7 > thresholds(ii));
    passes_all (ii) = sum ( (1-datasummary(:,7)) > thresholds(ii));
    %passes_LOOCV8 (ii) = sum ( dice_LOOCV8 > thresholds(ii));
end

%passes_iter_AUC = sum (passes_iter) ./ (10001 * stats_mu_iter.n) ;  % The AUC is actually the same as the mean
%passes_LOOCV7_AUC = sum (passes_LOOCV7) ./ (10001 * stats_mu7.n);
% passes_LOOCV7_AUC_Low = sum (passes_LOOCV_Low7) ./ (10001 * stats_muLow7.n);
% passes_LOOCV7_AUC_High = sum (passes_LOOCV_High7) ./ (10001 * stats_muHigh7.n);
%passes_LOOCV8_AUC = sum (passes_LOOCV8) ./ (10001 * stats_mu8.n);
%figure(1); plot (thresholds, passes_LOOCV8,'LineWidth',5);
%figure; plot (thresholds, passes_LOOCV7);
figure; plot (thresholds, passes_LOOCV25);
% figure; plot (thresholds, passes_LOOCV_Low7);
% figure; plot (thresholds, passes_LOOCV_High7);
figure; plot (thresholds, passes_all);
%figure(3); hist (mu_eff8);
%figure(4); hist (mu_eff7);

stats_dice25
