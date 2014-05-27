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
close all
hold off
load data.mat

% Load the optimized data
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
cell_data = csvimport('alt_datasummary.txt');
headers = cell_data(1,1:3);
data.opt = cell2mat(cell_data(2:end,:));
data.opt(:,2) = data.opt(:,2)*(6000-100)+100; % Convert from DAKOTA's scaled to absolute mu_eff


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
threshold7 = max_dice;
index_7= max_dice;
for ii = 1:num_runs
    for jj = 1:num_therms
        [max_dice(ii,jj),index_dice(ii,jj)] = max(data.dice(:,jj,ii)); % Answers Q1
        %[threshold7(ii,jj),index_7(ii,jjMEWCD4)] = find(data.dice(:,jj,ii)>0.7);
    end
    [min_obj_fxn(ii),index_obj_fxn(ii)] = min(data.data(:,5,ii));
    diff_dice(ii) = data.data(index_obj_fxn(ii),1,ii)-data.data(index_dice(ii),1,ii); % Answers Q2
    diff_opt_to_mu_eff(ii) = data.opt(ii,2) - data.data(index_obj_fxn(ii),1,ii); % Answers Q3a
    diff_opt_to_dice(ii) = data.opt(ii,2) - data.data(index_dice(ii),1,ii); % Answers Q3b
end

clear ii jj
dice_pass_count = zeros(num_therms,1);
dice_pass_7 = zeros(num_therms,num_runs);
index_pass_7 = dice_pass_7;
dice_pass_8 = dice_pass_7;
index_pass_8 = dice_pass_7;
dice_pass_83 = dice_pass_7;
index_pass_83=dice_pass_7;
for ii = 1:num_therms
    dice_pass_count(ii) = sum( max(data.dice(:,ii,:) > 0.7));
    for jj = 1:num_runs
        [dice_pass_7(ii,jj), index_pass_7(ii,jj)] = max(data.dice(:,ii,jj) >= 0.7);
        [dice_pass_8(ii,jj), index_pass_8(ii,jj)] = max(data.dice(:,ii,jj) >= 0.8);
        [dice_pass_83(ii,jj), index_pass_83(ii,jj)] = max(data.dice(:,ii,jj) >= 0.83);
    end
end
clear ii jj

% Find the Dice > 0.7 that also doesn't have indexes from the origin (ie
% near mu_eff = 0)
filter7 = index_pass_7 >=1000;  % index_pass_7 is for Dice >0.7
index_pass_7=index_pass_7.*filter7;
isotherm51_pass = (index_pass_7(1,:) >0).*index_pass_7(1,:);
isotherm57_pass = (index_pass_7(7,:) >0).*index_pass_7(7,:);
isotherm65_pass = (index_pass_7(15,:) >0).*index_pass_7(15,:);
mu_eff_pass51 = zeros (sum(isotherm51_pass > 0),3);
mu_eff_pass57 = zeros (sum(isotherm57_pass > 0),3);
mu_eff_pass65 = zeros (sum(isotherm65_pass > 0),3);
mu_eff_pass51(:,1) = isotherm51_pass(isotherm51_pass>0);
mu_eff_pass51(:,2) = find(isotherm51_pass>0);
mu_eff_pass51(:,3) = data.data(mu_eff_pass51(:,1),1,1);
mu_eff_pass57(:,1) = isotherm57_pass(isotherm57_pass>0);
mu_eff_pass57(:,2) = find(isotherm57_pass>0);
mu_eff_pass57(:,3) = data.data(mu_eff_pass57(:,1),1,1);
mu_eff_pass65(:,1) = isotherm65_pass(isotherm65_pass>0);
mu_eff_pass65(:,2) = find(isotherm65_pass>0);
mu_eff_pass65(:,3) = data.data(mu_eff_pass65(:,1),1,1);

filter8 = index_pass_8 >=1000;  % index_pass_7 is for Dice >0.7
index_pass_8=index_pass_8.*filter8;
isotherm51_pass8 = (index_pass_8(1,:) >0).*index_pass_8(1,:);
isotherm57_pass8 = (index_pass_8(7,:) >0).*index_pass_8(7,:);
isotherm65_pass8 = (index_pass_8(15,:) >0).*index_pass_8(15,:);
mu_eff_8pass51 = zeros (sum(isotherm51_pass8 > 0),3);
mu_eff_8pass57 = zeros (sum(isotherm57_pass8 > 0),3);
mu_eff_8pass65 = zeros (sum(isotherm65_pass8 > 0),3);
mu_eff_8pass51(:,1) = isotherm51_pass8(isotherm51_pass8>0);
mu_eff_8pass51(:,2) = find(isotherm51_pass8>0);
mu_eff_8pass51(:,3) = data.data(mu_eff_8pass51(:,1),1,1);
mu_eff_8pass57(:,1) = isotherm57_pass8(isotherm57_pass8>0);
mu_eff_8pass57(:,2) = find(isotherm57_pass8>0);
mu_eff_8pass57(:,3) = data.data(mu_eff_8pass57(:,1),1,1);
mu_eff_8pass65(:,1) = isotherm65_pass8(isotherm65_pass8>0);
mu_eff_8pass65(:,2) = find(isotherm65_pass8>0);
mu_eff_8pass65(:,3) = data.data(mu_eff_8pass65(:,1),1,1);

filter83 = index_pass_83 >=1000;  % index_pass_7 is for Dice >0.7
index_pass_83=index_pass_83.*filter83;
isotherm51_pass83 = (index_pass_83(1,:) >0).*index_pass_83(1,:);
isotherm57_pass83 = (index_pass_83(7,:) >0).*index_pass_83(7,:);
isotherm65_pass83 = (index_pass_83(15,:) >0).*index_pass_83(15,:);
mu_eff_83pass51 = zeros (sum(isotherm51_pass83 > 0),3);
mu_eff_83pass57 = zeros (sum(isotherm57_pass83 > 0),3);
mu_eff_83pass65 = zeros (sum(isotherm65_pass83 > 0),3);
mu_eff_83pass51(:,1) = isotherm51_pass83(isotherm51_pass83>0);
mu_eff_83pass51(:,2) = find(isotherm51_pass83>0);
mu_eff_83pass51(:,3) = data.data(mu_eff_83pass51(:,1),1,1);
mu_eff_83pass57(:,1) = isotherm57_pass83(isotherm57_pass83>0);
mu_eff_83pass57(:,2) = find(isotherm57_pass83>0);
mu_eff_83pass57(:,3) = data.data(mu_eff_83pass57(:,1),1,1);
mu_eff_83pass65(:,1) = isotherm65_pass83(isotherm65_pass83>0);
mu_eff_83pass65(:,2) = find(isotherm65_pass83>0);
mu_eff_83pass65(:,3) = data.data(mu_eff_83pass65(:,1),1,1);

% Show figures
% figure; hist ( mu_eff_pass51(:,3));
 figure; hist ( mu_eff_pass57(:,3));
% figure; hist ( mu_eff_pass65(:,3));

% figure; hist ( mu_eff_8pass51(:,3));
figure; hist ( mu_eff_8pass57(:,3));
% figure; hist ( mu_eff_8pass65(:,3));

% [ hh7_57 dice_7values57 hh7_57_alt dice_7values57_alt] = LOOCV_script_dice_input (mu_eff_pass57);
% [ hh8_57 dice_8values57 hh8_57_alt dice_8values57_alt] = LOOCV_script_dice_input (mu_eff_8pass57);
% [ hh83_57 dice_83values57 hh83_57_alt dice_83values57_alt] = LOOCV_script_dice_input (mu_eff_83pass57);
% figure; hist ( mu_eff_83pass51(:,3));
% figure; hist ( mu_eff_83pass57(:,3));
% figure; hist ( mu_eff_83pass65(:,3));


% Optimized mu_eff histograms
% figure; hist ( data.opt(:,2) );  % put in presentation
% stats.data_opt.mean = mean( data.opt(:,2));
% stats.data_opt.median = median( data.opt(:,2));
% stats.data_opt.std = std( data.opt(:,2));
% stats.data_opt.skew = skewness( data.opt(:,2)); % Positive is for long right tail (positivie tail); and vice versa; alternative postive skew is for bulk of distribution on left.
% stats.data_opt.kurt = kurtosis( data.opt(:,2)); % Positive is extra peaked-ness compared to Normal dist. and vice versa (ie neg means broader than Normal dist)
% temporary = data.opt(:,2);
% temporary = temporary(temporary > 11);
% stats.temp.mean = mean ( temporary );
% stats.temp.median = median ( temporary);
% stats.temp.std = std( temporary );
% stats.temp.skew = skewness ( temporary );
% stats.temp.kurt = kurtosis ( temporary );
% figure; hist (temporary); % put in presentation
% clear temporary

% Dice from optimized runs.
% cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
% aaa = load('opt_dice');
% hist_x_values = linspace(0.05,0.95,10);
% figure; hist (aaa.total(:,6),hist_x_values); % put in presentation
% stats.dice_opt.mean = mean( aaa.total(:,6));
% stats.dice_opt.median = median( aaa.total(:,6));
% stats.dice_opt.std = std ( aaa.total(:,6));
% stats.dice_opt.skew = skewness( aaa.total(:,6));
% stats.dice_opt.kurt = kurtosis( aaa.total(:,6));
% temporary = aaa.total(:,6) ; % Eliminate the runs where it optimized to mu_eff = 10;
% list = find ( aaa.total(:,1) > 11);
% temporary = temporary (list);
% figure; hist (temporary,hist_x_values); % put in presentation
% stats.dice_temp.mean = mean (temporary);
% stats.dice_temp.median = median( temporary );
% stats.dice_temp.std = std( temporary);
% stats.dice_temp.skew = skewness( temporary);
% stats.dice_temp.kurt = kurtosis( temporary);
% clear temporary hist_x_values

% mu_eff and dice from brute force sensitivity studies
% figure; plot(data.data(:,1,1),data.data(:,5,1));
% xlim([0,6000]);
% ylim([0,4*10^5]);
% figure; plot(data.data(1:2501,1,1),data.dice(1:2501,7,1));
% xlim([0,5000]);
% ylim([0.1,0.9]);
% figure; plot(data.data(1:2501,1,1),data.dice(1:2501,7,1),'b-');
% xlim([0,5000]);
% ylim([0.1,0.9]);
% hold on
% plot(data.data(1:2501,1,1),data.dice(1:2501,1,1),'r-');
% plot(data.data(1:2501,1,1),data.dice(1:2501,15,1),'g-');
% hold off
figure; hist (data.data(index_obj_fxn,1,1));

% Max Dice Histograms
% figure; hist (data.data(index_dice(:,1),1,1));
% figure; hist (data.data(index_dice(:,7),1,1));
% figure; hist (data.data(index_dice(:,15),1,1));

% % Drop the bad row (23)
%clean_index_dice = index_dice(1:22,:);
%clean_index_dice(23:24,:) = index_dice(24:25,:);
%figure; hist (data.data(clean_index_dice(:,1),1,1));
%figure; hist (data.data(clean_index_dice(:,7),1,1));
%figure; hist (data.data(clean_index_dice(:,15),1,1));


% How did the optimization perform? (diff of opt and brute force mu_eff)
% [~,index] = min( data.data(:,5,:));
% index = squeeze (index);
% brute.mu_eff = data.data(index,1,:);
% brute.mu_eff = brute.mu_eff(:,1,1);
% diff.mu_eff = aaa.total (:,1) - brute.mu_eff;
% hist_x_values = linspace(-1,1,11);
% figure; hist (diff.mu_eff(list),hist_x_values); % shows hist of diff, but without the 2 extremely wrong ones

% What is relationship between the obj_fxn and Dice? (diff of
% mu_eff@min_obj_fxn and mu_eff@maxDice)
%diff.diceVobjfxn = data.dice(:,7,:) - data.data(


% Dice from brute force sensitivity studies.
figure; plot (data.temp_headers,dice_pass_count,'b--o','MarkerSize',15);
