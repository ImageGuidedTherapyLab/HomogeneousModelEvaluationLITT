% Script for making plots for June 2014 committee meeting

% Make sure the opttype and Matlab_flag are set correctly
hist_plotting;

% Gathers all of the data into one nice group of variables
[L2norm, dice, tmap_model, MRTI_crop,mat_struct] = check_opt_GoodBadUgly( Study_paths, 'bestfit1', datasummary(:,3));

% For Aim 2, I want to show Good, Bad, and Ugly examples of the
% optimization. This is before the DSC > 0.8 threshold.

% This runs only Study0035/0530 ( Good example)
% Displays the model and MRTI temperatures and 57 C label mats. Then it
% also shows the physical dimension of the ablation.
[one_dice_Study35_530, model_Study35_530, MRTI_Study35_530,intersect_Study35_530] = single_iter_Dice(tmap_model{1},MRTI_crop{1});
(mat_struct(1).file.inputdatavars.voi(2)-mat_struct(1).file.inputdatavars.voi(1)+1) * mat_struct(1).file.inputdatavars.spacing(1) % This returns zero '0' for some frustrating reason
(mat_struct(1).file.inputdatavars.voi(4)-mat_struct(1).file.inputdatavars.voi(3)+1) * mat_struct(1).file.inputdatavars.spacing(2)

% Now I want to find a Bad one. It'll probably be between 0.1 and 0.6
aaa=(dice_raw < 0.6).*(dice_raw > 0.1); % Find satisfactory datasets
bbb=find(aaa==0);
ccc=dice_raw; % Make dummy list
ccc(bbb) = []; % Toss unwanted data.
Bad_paths = Study_paths;
Bad_paths(bbb,:) = [];
Bad_tmap = tmap_model;
Bad_tmap(bbb) = [];
Bad_MRTI = MRTI_crop;
Bad_MRTI(bbb) = [];
Bad_struct = mat_struct;
Bad_struct(bbb) = [];

%[one_dice_Study25_436, model_Study25_436,MRTI_Study24_436]=single_iter_Dice(Bad_tmap{1},Bad_MRTI{1});%Poor registration

[one_dice_Study28_468, model_Study28_468, MRTI_Study28_468] = single_iter_Dice(Bad_tmap{3},Bad_MRTI{3});
(Bad_struct(3).file.inputdatavars.voi(2)-Bad_struct(3).file.inputdatavars.voi(1)+1) * Bad_struct(3).file.inputdatavars.spacing(1) % This returns zero '0' for some frustrating reason
(Bad_struct(3).file.inputdatavars.voi(4)-Bad_struct(3).file.inputdatavars.voi(3)+1) * Bad_struct(3).file.inputdatavars.spacing(2)


[one_dice_Study28_471, model_Study28_471, MRTI_Study28_471] = single_iter_Dice(Bad_tmap{4},Bad_MRTI{4});
(Bad_struct(4).file.inputdatavars.voi(2)-Bad_struct(4).file.inputdatavars.voi(1)+1) * Bad_struct(4).file.inputdatavars.spacing(1) % This returns zero '0' for some frustrating reason
(Bad_struct(4).file.inputdatavars.voi(4)-Bad_struct(4).file.inputdatavars.voi(3)+1) * Bad_struct(4).file.inputdatavars.spacing(2)

clear aaa bbb ccc
% This following section shows that the optimization finds the global
% maximum by comparing to a direct search algorithm
aaa = load ( '/FUS4/data2/sjfahrenholtz/MATLAB/Tests/June_committee_direct_search_large.mat');
figure(9); plot( aaa.total_path{1,2}(1:5000,1), aaa.total_path{1,3}(1:5000,7));
[dice_max , index_max] = max( aaa.total_path{1,3}(:,7)); % Compare the direct search's max to the pattern search
multi_maxes = find( aaa.total_path{1,3}(:,7) == dice_max); % There are multiple maxima for the direct search

% This explores why there is a shoulder on the DSC

Brute_sensitivity22;

% This section shows the extrema of the Dice from the Direct search
figure(9); plot( aaa.total_path{1,2}(1:5000,1), aaa.total_path{1,3}(1:5000,1), aaa.total_path{1,2}(1:5000,1), aaa.total_path{1,3}(1:5000,7), aaa.total_path{1,2}(1:5000,1), aaa.total_path{1,3}(1:5000,15));
% Check where the max are at
[dice_max51 , index_max51] = max( aaa.total_path{1,3}(:,1));
[dice_max57 , index_max57] = max( aaa.total_path{1,3}(:,7));
[dice_max65 , index_max65] = max( aaa.total_path{1,3}(:,15));


% This section shows the plausible isotherm choices
figure(9); plot( aaa.total_path{1,2}(1:5000,1), aaa.total_path{1,3}(1:5000,6), aaa.total_path{1,2}(1:5000,1), aaa.total_path{1,3}(1:5000,10));
[dice_max56 , index_max56] = max( aaa.total_path{1,3}(:,6));
[dice_max60 , index_max60] = max( aaa.total_path{1,3}(:,10));

% This section summarizes the overall performance of the pattern search
% optimization and is the last section for Aim 2.
thresholds = linspace ( 0.0, 0.9, 10 );
passing_number = zeros(10,1);
delta_number = zeros(10,1);

for ii = 1:10
    passing_number (ii) = sum ( (1-datasummary (:,7) ) >= thresholds (ii) );
    if  (ii ~= 1)
        delta_number(ii) = passing_number(ii-1) - passing_number(ii);
    end
end
figure(10),plot(thresholds,passing_number, 'LineWidth', 5)


% This section begins Aim 3
figure(3),hist(mu_eff8,15)

% Show the LOOCV performance of Study0030/0497 which is index 3
aaa=(dice_raw >= 0.8); % Find satisfactory datasets
bbb=find(aaa==0);
ccc=dice_raw; % Make dummy list
ccc(bbb) = []; % Toss unwanted data.
LOOCV_paths = Study_paths;
LOOCV_paths(bbb,:) = [];
LOOCV_tmap = tmap_model;
LOOCV_tmap(bbb) = [];
LOOCV_MRTI = MRTI_crop;
LOOCV_MRTI(bbb) = [];
LOOCV_struct = mat_struct;
LOOCV_struct(bbb) = [];


mu_eff_LOOCV = mu_eff8; % Make a copy of both the mu_eff values and the paths

[L2norm, dice, tmap_LOOCV, MRTI_LOOCV,mat_struct_LOOCV] = check_opt_LOOCV( LOOCV_paths, 'bestfit1', mu_eff_LOOCV );
[ ~, model_deg_threshold, MRTI_deg_threshold, intersection ] = single_iter_Dice (tmap_LOOCV{3}, MRTI_LOOCV{3});
(LOOCV_struct(3).file.inputdatavars.voi(2)-LOOCV_struct(3).file.inputdatavars.voi(1)+1) * LOOCV_struct(3).file.inputdatavars.spacing(1) % This returns zero '0' for some frustrating reason
(LOOCV_struct(3).file.inputdatavars.voi(4)-LOOCV_struct(3).file.inputdatavars.voi(3)+1) * LOOCV_struct(3).file.inputdatavars.spacing(2)

% Study0035/0530
[ ~, model_deg_threshold, MRTI_deg_threshold, intersection ] = single_iter_Dice (tmap_LOOCV{1}, MRTI_LOOCV{1});
(LOOCV_struct(1).file.inputdatavars.voi(2)-LOOCV_struct(1).file.inputdatavars.voi(1)+1) * LOOCV_struct(1).file.inputdatavars.spacing(1) % This returns zero '0' for some frustrating reason
(LOOCV_struct(1).file.inputdatavars.voi(4)-LOOCV_struct(1).file.inputdatavars.voi(3)+1) * LOOCV_struct(1).file.inputdatavars.spacing(2)

% Study0026/0450
[ ~, model_deg_threshold, MRTI_deg_threshold, intersection ] = single_iter_Dice (tmap_LOOCV{9}, MRTI_LOOCV{9});
(LOOCV_struct(9).file.inputdatavars.voi(2)-LOOCV_struct(9).file.inputdatavars.voi(1)+1) * LOOCV_struct(9).file.inputdatavars.spacing(1) % This returns zero '0' for some frustrating reason
(LOOCV_struct(9).file.inputdatavars.voi(4)-LOOCV_struct(9).file.inputdatavars.voi(3)+1) * LOOCV_struct(9).file.inputdatavars.spacing(2)

% Study0026/0414
[ ~, model_deg_threshold, MRTI_deg_threshold, intersection ] = single_iter_Dice (tmap_LOOCV{10}, MRTI_LOOCV{10});
(LOOCV_struct(10).file.inputdatavars.voi(2)-LOOCV_struct(10).file.inputdatavars.voi(1)+1) * LOOCV_struct(10).file.inputdatavars.spacing(1) % This returns zero '0' for some frustrating reason
(LOOCV_struct(10).file.inputdatavars.voi(4)-LOOCV_struct(10).file.inputdatavars.voi(3)+1) * LOOCV_struct(10).file.inputdatavars.spacing(2)

%DF vs SJF DVH

DF_results = load('/FUS4/data2/sjfahrenholtz/MATLAB/Tests/CombineDVHPlot.mat');
figure(2); plot (thresholds, passes_LOOCV8, thresholds, DF_results.passes_LOOCV8,'LineWidth',5);

DF_results.dice_LOOCV8_stats = Descriptive_statistics(DF_results.dice_LOOCV8);
DF_results.LOOCV8_passes7 = sum( DF_results.dice_LOOCV8 > 0.7);
dice_LOOCV8_stats = Descriptive_statistics(dice_LOOCV8);
LOOCV8_passes7 = sum( dice_LOOCV8 > 0.7);

[hh_comparison, p_comparison, CI_comparison, stats_comparison] =ttest2 ( dice_LOOCV8, DF_results.dice_LOOCV8)
