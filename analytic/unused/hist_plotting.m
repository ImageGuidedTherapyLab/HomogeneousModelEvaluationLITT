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


% Identify the studies to be examined.
% Study_paths {1,1} = 'Study0035';
% Study_paths {1,2} = '0530';
% Study_paths {2,1} = 'Study0023';
% Study_paths {2,2} = '0430';
% Study_paths {3,1} = 'Study0023';
% Study_paths {3,2} = '0428';
% Study_paths {4,1} = 'Study0030';
% Study_paths {4,2} = '0495';
% Study_paths {5,1} = 'Study0030';
% Study_paths {5,2} = '0497';
% Study_paths {6,1} = 'Study0030';
% Study_paths {6,2} = '0488';
% Study_paths {7,1} = 'Study0030';
% Study_paths {7,2} = '0491';
% Study_paths {8,1} = 'Study0030';
% Study_paths {8,2} = '0496';
% Study_paths {9,1} = 'Study0030';
% Study_paths {9,2} = '0490';
% Study_paths {10,1} = 'Study0017';
% Study_paths {10,2} = '0378';
% Study_paths {11,1} = 'Study0018';
% Study_paths {11,2} = '0402';
% Study_paths {12,1} = 'Study0018';
% Study_paths {12,2} = '0389';
% Study_paths {13,1} = 'Study0018';
% Study_paths {13,2} = '0385';
% Study_paths {14,1} = 'Study0029';
% Study_paths {14,2} = '0476';
% Study_paths {15,1} = 'Study0029';
% Study_paths {15,2} = '0477';
% 
% Study_paths {8,1} = 'Study0025';
% Study_paths {8,2} = '0438';
% Study_paths {9,1} = 'Study0025';
% Study_paths {9,2} = '0435';
% Study_paths {10,1} = 'Study0025';
% Study_paths {10,2} = '0440';
% Study_paths {11,1} = 'Study0025';
% Study_paths {11,2} = '0436';
% Study_paths {12,1} = 'Study0028';
% Study_paths {12,2} = '0466';
% Study_paths {13,1} = 'Study0028';
% Study_paths {13,2} = '0468';
% Study_paths {14,1} = 'Study0028';
% Study_paths {14,2} = '0471';
% Study_paths {15,1} = 'Study0026';
% Study_paths {15,2} = '0447';
% Study_paths {16,1} = 'Study0026';
% Study_paths {16,2} = '0457';
% Study_paths {17,1} = 'Study0026';
% Study_paths {17,2} = '0455';
% Study_paths {18,1} = 'Study0026';
% Study_paths {18,2} = '0453';
% Study_paths {19,1} = 'Study0026';
% Study_paths {19,2} = '0450';
% Study_paths {20,1} = 'Study0026';
% Study_paths {20,2} = '0451';
% Study_paths {21,1} = 'Study0022';
% Study_paths {21,2} = '0418';
% Study_paths {22,1} = 'Study0022';
% Study_paths {22,2} = '0417';
% Study_paths {23,1} = 'Study0021';
% Study_paths {23,2} = '0409';
% Study_paths {24,1} = 'Study0021';
% Study_paths {24,2} = '0414';
% Study_paths {25,1} = 'Study0021';
% Study_paths {25,2} = '0415';


% read  best_fit optimization data and store mu_eff and alpha
%datasummary = dlmread('ex_datasummary.txt',',',1,0);
datasummary = dlmread('datasummary.txt',',',1,0);
datasummary(any(isnan(datasummary), 2), :) = [];
num_studies = size(datasummary,1);

Study_paths = cell (num_studies,2);
for ii = 1:num_studies
    
    Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end

% best_iter = datasummary(:,3);
% mu_eff_data = datasummary(:,4);
% alpha_data = datasummary(:,5);
% dice_opt = datasummary(:,7);
% L2norm_opt = datasummary(:,8);




% num_datasummary = size(datasummary,1);
% 
% matching_num = zeros(1,num_studies);
% mu_eff_index = zeros(1,num_studies);
% mu_eff_opt   = zeros(1,num_studies);
% copy_summary = zeros(num_studies,8);
% for ii = 1:num_studies
%     matching_num(ii) = str2num(Study_paths{(ii),2});
%     mu_eff_index(ii) = find( datasummary(:,1) == matching_num(ii));
%     copy_summary(ii,:) = datasummary( mu_eff_index(ii) , : );
%     %mu_eff_opt(ii) = dlm_data(mu_eff_index(ii),2);
% end
% clear ii
% 
% copy_summary( isnan( copy_summary ) ) = 0;
% best_iter   = copy_summary(:,2);
% mu_eff_data = copy_summary(:,3);
% alpha_data  = copy_summary(:,4);

% TODO testing for now

% matching_num = zeros(1,num_datasummary);
% mu_eff_index = matching_num;
% mu_eff_opt   = matching_num;
% for ii = 1:num_studies
%     matching_num(ii) = str2num(Study_paths{(ii),2});
%     mu_eff_index(ii) = find( datasummary(:,1) == matching_num(ii));
%     mu_eff_opt(ii) = data_summary(mu_eff_index(ii));
% end
% clear ii



% for ii = 1:num_studies
%     if  isnan(dlm_data(ii,3)) == 1
%         dlm_data(ii,3) =  1;
%     end
% end
% clear ii

% matching_num = zeros(1,num_studies);
% mu_eff_index = zeros(1,num_studies);
% mu_eff_opt   = zeros(1,num_studies);
% for ii = 1:num_studies
%     matching_num(ii) = str2num(Study_paths{(ii),2});
%     mu_eff_index(ii) = find( dlm_data(:,1) == matching_num(ii));
%     mu_eff_opt(ii) = dlm_data(mu_eff_index(ii),2);
% end
% clear ii

%mu_eff_opt22 = mu_eff_opt*(6000-100)+100; %This line converts the normalized value into absolute. It's very import to get the converion correct

% This script finds the best mu_eff for the different studies.
opttype = 'bestfit' ;
dice_raw = datasummary(:,7);
%opttype = 'bestfit1' ;
%dice_raw = 1 - datasummary(:,7);

toss_index7 = find(dice_raw<0.7);
toss_index8 = find(dice_raw<0.8);

mu_eff7 = datasummary(:,4);
mu_eff7 (toss_index7) = [];
mu_eff8 = datasummary(:,4);
mu_eff8 (toss_index8) = [];

alpha7 = datasummary(:,5);
alpha7 (toss_index7) = [];
alpha8 = datasummary(:,5);
alpha8 (toss_index8) = [];

best_iter7 = datasummary(:,3);
best_iter7 (toss_index7) = [];
best_iter8 = datasummary(:,3);
best_iter8 (toss_index8) = [];

dice7 = dice_raw;
dice7(toss_index7) = [];
dice8 = dice_raw;
dice8(toss_index8) = [];

stats_mu_raw = Descriptive_statistics( datasummary(:,4) );
stats_alpha_raw = Descriptive_statistics( datasummary(:, 5) );
stats_mu7 = Descriptive_statistics(mu_eff7);
stats_alpha7 = Descriptive_statistics(alpha7);
stats_mu8 = Descriptive_statistics(mu_eff8);
stats_alpha8 = Descriptive_statistics(alpha8);

% figure; hist(mu_eff_opt22);
% figure; hist(mu_eff7);
% figure; hist(mu_eff8);

% Remove Study_paths indices
temp_paths = Study_paths;
temp_paths(toss_index7,:) = [];
Study_paths7 = temp_paths;
temp_paths = Study_paths;
temp_paths(toss_index8,:) = [];
Study_paths8 = temp_paths;

%[ hh_raw, dice_values_LOOCV ] = LOOCV_t_test_DiceTemp ( Study_paths, datasummary(:,4) ,datasummary(:,5) , datasummary(:,3), opttype );
%[ hh7, dice_LOOCV7 ] = LOOCV_t_test_DiceTemp ( Study_paths7, mu_eff7, alpha7, best_iter7, opttype );
[ hh8, dice_LOOCV8 ] = LOOCV_t_test_DiceTemp ( Study_paths8, mu_eff8, alpha8, best_iter8, opttype );

%hh8.ptest = .5;

% mu_eff_iter = mu_eff7;
% paths_iter = Study_paths7;
% alpha_iter = alpha7;
% best_iter_iter = best_iter7;
% hh_iter.ptest = hh7.ptest;
% stats_mu_iter = stats_mu7;
% iteration_tracker = stats_mu_iter.n;

% mu_eff_iter = mu_eff8;
% paths_iter = Study_paths8;
% alpha_iter = alpha8;
% best_iter_iter = best_iter8;
% hh_iter.ptest = hh8.ptest;
% stats_mu_iter = stats_mu8;
% iteration_tracker = stats_mu_iter.n;
% 
% while hh_iter.ptest > 0.05 && stats_mu_iter.n > 3
%     
%     disp( iteration_tracker );
%     disp( hh_iter.ptest );
%     residual8 = abs( mu_eff_iter - stats_mu_iter.mean ); % Find the largest residual
%     [~,toss_index] = max( residual8 );
%     mu_eff_iter(toss_index) = [];
%     alpha_iter (toss_index) = [];
%     paths_iter (toss_index,:) = [];
%     best_iter_iter (toss_index) = [];
%     stats_mu_iter = Descriptive_statistics( mu_eff_iter);
%     iteration_tracker = stats_mu_iter;
%     
%     [ hh_iter, dice_values_iter ] = LOOCV_t_test_DiceTemp ( paths_iter, mu_eff_iter, alpha_iter, best_iter_iter, opttype );
%     
% end

% if hh_iter.ptest < 0.05 && stats_mu_iter.n > 3
%     stats_mu_iter = Descriptive_statistics( mu_eff_iter);
%     stats_alpha_iter = Descriptive_statistics( alpha_iter );
% end

thresholds = linspace ( 0.0, 1, 10001);
passes_LOOCV8 = zeros (10001,1);
%passes_LOOCV7 = passes_LOOCV8;
%passes_iter = passes_LOOCV8;
for ii = 1:10001
    %passes_iter   (ii) = sum ( dice_values_iter > thresholds(ii));
    %passes_LOOCV7 (ii) = sum ( dice_LOOCV7 > thresholds(ii));
    passes_LOOCV8 (ii) = sum ( dice_LOOCV8 > thresholds(ii));
end

%passes_iter_AUC = sum (passes_iter) ./ (10001 * stats_mu_iter.n) ;  % The AUC is actually the same as the mean
%passes_LOOCV7_AUC = sum (passes_LOOCV7) ./ (10001 * stats_mu7.n);
passes_LOOCV8_AUC = sum (passes_LOOCV8) ./ (10001 * stats_mu8.n);
figure(1); plot (thresholds, passes_LOOCV8);
%figure(2); plot (thresholds, passes_LOOCV7);
figure(3); hist (mu_eff8);
%figure(4); hist (mu_eff7);
