% This script is meant to see how many datasets have solutions that are
% good. The following questions need to be answered:
% 1. Does a good solution exist? This will be done by doing max(dice) >0.7
% 2. How far (in mu_eff value) is the max Dice from the min obj fxn?
%              I.e. (mu_eff@maxDice - mu_eff@minObjFxn)
% 3. Where does the optimizer go? Ie (mu_eff@Opt - mu_eff@minObjFxn) and
%           (mu_eff@Opt - mu_eff@maxDice)
% 4. How does the maxDice depend on isotherm choice (temperature).

cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
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
load data.mat

% Load the optimized data
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
cell_data = csvimport('alt_datasummary.txt');
headers = cell_data(1,1:3);
data.opt = cell2mat(cell_data(2:end,:));
data.opt(:,2) = data.opt(:,2)*(15000-10)+10; % Convert from DAKOTA's scaled to absolute mu_eff


% Question 1
num_runs = size(data.labels,1);
num_therms = size(data.dice,2);
max_dice = zeros (num_runs,1);
index_dice = max_dice;
min_obj_fxn = max_dice;
index_obj_fxn = max_dice;
diff_dice = max_dice;
diff_opt_to_mu_eff = max_dice;
diff_opt_to_dice = max_dice;
for ii = 1:num_runs
    for jj = 1:num_therms
        [max_dice(ii,jj),index_dice(ii,jj)] = max(data.dice(:,jj,ii)); % Answers Q1
    end
    [min_obj_fxn(ii),index_obj_fxn(ii)] = min(data.data(:,5,ii));
    diff_dice(ii) = data.data(index_obj_fxn(ii),1,ii)-data.data(index_dice(ii),1,ii); % Answers Q2
    diff_opt_to_mu_eff(ii) = data.opt(ii,2) - data.data(index_obj_fxn(ii),1,ii); % Answers Q3a
    diff_opt_to_dice(ii) = data.opt(ii,2) - data.data(index_dice(ii),1,ii); % Answers Q3b
end

clear ii jj
dice_pass_count = zeros(num_therms,1);
for ii = 1:num_therms
    dice_pass_count(ii) = sum( max(data.dice(:,ii,:) > 0.7));
end
clear ii
% Optimized mu_eff histograms
figure(1); hist ( data.opt(:,2) );  % put in presentation
stats.data_opt.mean = mean( data.opt(:,2));
stats.data_opt.median = median( data.opt(:,2));
stats.data_opt.std = std( data.opt(:,2));
stats.data_opt.skew = skewness( data.opt(:,2)); % Positive is for long right tail (positivie tail); and vice versa; alternative postive skew is for bulk of distribution on left.
stats.data_opt.kurt = kurtosis( data.opt(:,2)); % Positive is extra peaked-ness compared to Normal dist. and vice versa (ie neg means broader than Normal dist)
temporary = data.opt(:,2);
temporary = temporary(temporary > 11);
stats.temp.mean = mean ( temporary );
stats.temp.median = median ( temporary);
stats.temp.std = std( temporary );
stats.temp.skew = skewness ( temporary );
stats.temp.kurt = kurtosis ( temporary );
figure(2); hist (temporary); % put in presentation
clear temporary
% Dice from optimized runs.
cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
aaa = load('opt_dice');
hist_x_values = linspace(0.05,0.95,10);
figure (3); hist (aaa.total(:,6),hist_x_values); % put in presentation
stats.dice_opt.mean = mean( aaa.total(:,6));
stats.dice_opt.median = median( aaa.total(:,6));
stats.dice_opt.std = std ( aaa.total(:,6));
stats.dice_opt.skew = skewness( aaa.total(:,6));
stats.dice_opt.kurt = kurtosis( aaa.total(:,6));
temporary = aaa.total(:,6) ; % Eliminate the runs where it optimized to mu_eff = 10;
list = find ( aaa.total(:,1) > 11);
temporary = temporary (list);
figure(4); hist (temporary,hist_x_values); % put in presentation
stats.dice_temp.mean = mean (temporary);
stats.dice_temp.median = median( temporary );
stats.dice_temp.std = std( temporary);
stats.dice_temp.skew = skewness( temporary);
stats.dice_temp.kurt = kurtosis( temporary);
clear temporary hist_x_values

% mu_eff and dice from brute force sensitivity studies
figure(5); plot(data.data(:,1),data.data(:,5,1));
figure(6); plot(data.data(:,1),data.dice(:,7,1));

% How did the optimization perform? (diff of opt and brute force mu_eff)
[~,index] = min( data.data(:,5,:));
index = squeeze (index);
brute.mu_eff = data.data(index,1,:);
brute.mu_eff = brute.mu_eff(:,1,1);
diff.mu_eff = aaa.total (:,1) - brute.mu_eff;
hist_x_values = linspace(-1,1,11);
figure(7); hist (diff.mu_eff(list),hist_x_values); % shows hist of diff, but without the 2 extremely wrong ones

% What is relationship between the obj_fxn and Dice? (diff of
% mu_eff@min_obj_fxn and mu_eff@maxDice)
diff.diceVobjfxn = data.dice(:,7,:) - data.data(


% Dice from brute force sensitivity studies.
figure(11); plot (data.temp_headers,dice_pass_count,'b--o','MarkerSize',15);
