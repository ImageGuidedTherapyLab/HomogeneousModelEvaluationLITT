function [stats, LOOCV_stats] = master_LOOCV ( opttype, data_filename, dice_thresholds, mu_thresholds, Matlab_flag);

datasummary = dlmread(data_filename,',',1,0);
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

% variable initialization

if isempty(dice_thresholds) ==1
    
    stats.dice.labels = cell(1,1);
    stats.dice.labels = 'opt';
    stats.dice.values = cell(1,1);
    length_dice_thresholds = 0;
else

    length_dice_thresholds = length(dice_thresholds);
    toss_dice = cell(length_dice_thresholds,1);
    for jj = 1:length_dice_thresholds
        
        toss_dice{jj}= find(dice_thresholds(jj) > datasummary(:,7));
        
    end
    
end

if isempty(mu_thresholds) ==1
    
    stats.mu_eff.labels = cell(1,1);
    stats.mu_eff.labels = 'opt';
    stats.mu_eff.values = cell(1,1);
    length_mu_groups = 0;
else
    
    length_mu_groups = length(mu_thresholds)+1;    
    toss_mu = cell(length_mu_groups,1);
    toss_mu{1} = find(datasummary(:,4) > mu_thresholds(1));
    toss_mu{end} = find(mu_thresholds(end) >= datasummary(:,4));
    
    for ii=2:(length_mu_groups-1)
        
        toss_mu{ii} = find(  (mu_thresholds(ii-1) >= datasummary(:,4)) + (datasummary(:,4) > mu_thresholds(ii))  );
        
    end
end

clear ii jj

% Permute the toss variables.
total_toss = cell(length_mu_groups, length_dice_thresholds);
mu_eff = total_toss;
alpha = total_toss;
best_iter = total_toss;
dice = total_toss;
LOOCV_dice = total_toss;
stats.dice.labels = total_toss;
stats.dice.values = total_toss;
stats.mu_eff.labels = total_toss;
stats.mu_eff.values = total_toss;
LOOCV_stats.dice.labels = total_toss;
LOOCV_stats.dice.values = total_toss;
LOOCV_stats.dice.hh=total_toss;
kk=1;
for ii = 1:length_mu_groups
    for jj = 1:length_dice_thresholds
        
        total_toss{ii,jj} = unique([ toss_mu{ii}; toss_dice{jj}]);
        
        mu_eff{ii,jj} = datasummary(:,4);
        mu_eff{ii,jj} (total_toss{ii,jj})=[];
        
        alpha{ii,jj} = datasummary(:,5);
        alpha{ii,jj} (total_toss{ii,jj}) = [];
        
        best_iter{ii,jj} = datasummary(:,3);
        best_iter{ii,jj} (total_toss{ii,jj})=[];
        
        dice{ii,jj} = datasummary(:,7);
        dice{ii,jj}(total_toss{ii,jj}) = [];
        
        stats.dice.values{ii,jj} = Descriptive_statistics(dice{ii,jj});
        stats.mu_eff.values{ii,jj}   = Descriptive_statistics(mu_eff{ii,jj});
        
        if ii ==1
            stats.dice.labels{ii,jj} = strcat(['mu_eff <'], num2str(mu_thresholds(ii)), ['   dice ='], num2str(dice_thresholds(jj)) );
        elseif ii == length_mu_groups
            stats.dice.labels{ii,jj} = strcat(['mu_eff >'], num2str(mu_thresholds(ii-1)), ['   dice ='], num2str(dice_thresholds(jj)) );     
        else
            stats.dice.labels{ii,jj} = strcat(['mu_eff ='], num2str(mu_thresholds(ii-1)), [' to '], num2str(mu_thresholds(ii)), ['   dice ='], num2str(dice_thresholds(jj)) );
        end
        
        stats.mu_eff.labels{ii,jj} = stats.dice.labels{ii,jj};
        LOOCV_stats.dice.labels{ii,jj} = stats.dice.labels{ii,jj};
        
        temp_paths = Study_paths;
        temp_paths(total_toss{ii,jj},:) = [];
        
        disp(stats.dice.labels{ii,jj});
        disp(num2str(kk), [' of '], num2str(length_mu_groups .* length_dice_thresholds));
        [ LOOCV_stats.dice.hh{ii,jj}, LOOCV_dice{ii,jj}] = LOOCV_t_test_DiceTemp( temp_paths, mu_eff{ii,jj}, alpha{ii,jj}, best_iter{ii,jj}, opttype, Matlab_flag); 
        
        LOOCV_stats.dice.values{ii,jj} = Descriptive_statistics( LOOCV_dice{ii,jj});
    end
end





Matlab_flag = 1; % 0 means use FEM kernel; 1 means use MATLAB for kernel




% read  best_fit optimization data and store mu_eff and alpha
%datasummary = dlmread('ex_datasummary.txt',',',1,0);
datasummary = dlmread(data_filename,',',1,0);
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

dice_raw = datasummary(:,7);

toss_index7 = find(dice_raw<0.7);
%toss_index25= find(dice_raw<0.25);
toss_index3 = find(dice_raw<0.3);
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
% mu_eff25 = datasummary(:,4);
% mu_eff25 (toss_index25) = [];
mu_eff3 = datasummary(:,4);
mu_eff3 (toss_index3)=[];
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
% alpha25 = datasummary(:,5);
% alpha25 (toss_index25) = [];
alpha3 = datasummary(:,5);
alpha3 (toss_index3)=[];
% alpha8 = datasummary(:,5);
% alpha8 (toss_index8) = [];
% alphaLow7 = alpha7;
% alphaLow7(toss_indexLow7)=[];
% alphaHigh7=alpha7;
% alphaHigh7(toss_indexHigh7)=[];

best_iter7 = datasummary(:,3);
best_iter7 (toss_index7) = [];
% best_iter25 = datasummary(:,3);
% best_iter25 (toss_index25) = [];
best_iter3=datasummary(:,3);
best_iter3(toss_index3)=[];
% best_iter8 = datasummary(:,3);
% best_iter8 (toss_index8) = [];
% best_iterLow7 = best_iter7;
% best_iterLow7(toss_indexLow7)=[];
% best_iterHigh7=best_iter7;
% best_iterHigh7(toss_indexHigh7)=[];

dice7 = dice_raw;
dice7(toss_index7) = [];
% dice25 = dice_raw;
% dice25(toss_index25) = [];
dice3= dice_raw;
dice3(toss_index3)=[];
% dice8 = dice_raw;
% dice8(toss_index8) = [];
% diceLow7 = dice7;
% diceLow7(toss_indexLow7)=[];
% diceHigh7=dice7;
% diceHigh7(toss_indexHigh7)=[];

quality7 = quality;
quality7(toss_index7) = [];
% quality25 = quality;
% quality25(toss_index25) = [];
quality3=quality;
quality3(toss_index3)=[];

stats_mu_raw = Descriptive_statistics( datasummary(:,4) );
stats_alpha_raw = Descriptive_statistics( datasummary(:, 5) );
stats_mu7 = Descriptive_statistics(mu_eff7);
% stats_mu25 = Descriptive_statistics(mu_eff25);
stats_mu3 =  Descriptive_statistics(mu_eff3);
% stats_muLow7 = Descriptive_statistics(mu_effLow7);
% stats_muHigh7 = Descriptive_statistics(mu_effHigh7);
stats_alpha7 = Descriptive_statistics(alpha7);
% stats_alpha25 = Descriptive_statistics(alpha25);
stats_alpha3 = Descriptive_statistics(alpha3);
% stats_mu8 = Descriptive_statistics(mu_eff8);
% stats_alpha8 = Descriptive_statistics(alpha8);
stats_dice_raw = Descriptive_statistics(dice_raw);
stats_dice7 = Descriptive_statistics(dice7);
% stats_dice25 = Descriptive_statistics(dice25);
stats_dice3 = Descriptive_statistics(dice3);
% stats_diceLow7 = Descriptive_statistics(diceLow7);
% stats_diceHigh7 = Descriptive_statistics(diceHigh7);


figure; hist(datasummary(:,4));
figure; hist(mu_eff7);
% figure; hist(mu_eff25);
figure; hist(mu_eff3);
% figure; hist(mu_effLow7);
% figure; hist(mu_effHigh7);

% Remove Study_paths indices
temp_paths = Study_paths;
temp_paths(toss_index7,:) = [];
Study_paths7 = temp_paths;
% temp_paths = Study_paths;
% temp_paths(toss_index25,:) = [];
% Study_paths25 = temp_paths;
temp_paths = Study_paths;
temp_paths(toss_index3,:) = [];
Study_paths3 = temp_paths;
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
[ hh7, dice_LOOCV7 ] = LOOCV_t_test_DiceTemp ( Study_paths7, mu_eff7, alpha7, best_iter7, opttype, Matlab_flag );
% [ hh25, dice_LOOCV25 ] = LOOCV_t_test_DiceTemp ( Study_paths25, mu_eff25, alpha25, best_iter25, opttype, Matlab_flag );
[ hh3, dice_LOOCV3 ] = LOOCV_t_test_DiceTemp ( Study_paths3, mu_eff3, alpha3, best_iter3, opttype, Matlab_flag );
% [ hhLow7, dice_LOOCV_Low7 ] = LOOCV_t_test_DiceTemp ( Study_pathsLow7, mu_effLow7, alphaLow7, best_iterLow7, opttype, Matlab_flag );
% [ hhHigh7, dice_LOOCV_High7 ] = LOOCV_t_test_DiceTemp ( Study_pathsHigh7, mu_effHigh7, alphaHigh7, best_iterHigh7, opttype, Matlab_flag );
%[ hh8, dice_LOOCV8 ] = LOOCV_t_test_DiceTemp ( Study_paths8, mu_eff8, alpha8, best_iter8, opttype, Matlab_flag );
%dice_LOOCV8_stats = Descriptive_statistics(dice_LOOCV8);
dice_LOOCV7_stats = Descriptive_statistics(dice_LOOCV7);
% dice_LOOCV25_stats = Descriptive_statistics(dice_LOOCV25);
dice_LOOCV3_stats = Descriptive_statistics(dice_LOOCV3);
% dice_LOOCV_Low7stats = Descriptive_statistics(dice_LOOCV_Low7);
% dice_LOOCV_High7stats = Descriptive_statistics(dice_LOOCV_High7);


thresholds = linspace ( 0.0, 1, 10001);
%passes_LOOCV8 = zeros (10001,1);
passes_LOOCV7 = zeros (10001,1);
% passes_LOOCV25 = zeros (10001,1);
passes_LOOCV3 = zeros(10001,1);
% passes_LOOCV_Low7 = zeros(10001,1);
% passes_LOOCV_High7 = zeros(10001,1);
passes_all = zeros(10001,1);
%passes_iter = passes_LOOCV8;
for ii = 1:10001
    %passes_iter   (ii) = sum ( dice_values_iter > thresholds(ii));
    passes_LOOCV7 (ii) = sum ( dice_LOOCV7 > thresholds(ii));
%     passes_LOOCV25 (ii) = sum ( dice_LOOCV25 > thresholds(ii));
    passes_LOOCV3 (ii) = sum ( dice_LOOCV3 > thresholds(ii));
%     passes_LOOCV_Low7 (ii) = sum ( dice_LOOCV_Low7 > thresholds(ii));
%     passes_LOOCV_High7 (ii) = sum ( dice_LOOCV_High7 > thresholds(ii));
    passes_all (ii) = sum ( datasummary(:,7) > thresholds(ii));
    %passes_LOOCV8 (ii) = sum ( dice_LOOCV8 > thresholds(ii));
end

%passes_iter_AUC = sum (passes_iter) ./ (10001 * stats_mu_iter.n) ;  % The AUC is actually the same as the mean
passes_LOOCV7_AUC = sum (passes_LOOCV7) ./ (10001 * stats_mu7.n);
% passes_LOOCV7_AUC_Low = sum (passes_LOOCV_Low7) ./ (10001 * stats_muLow7.n);
% passes_LOOCV7_AUC_High = sum (passes_LOOCV_High7) ./ (10001 * stats_muHigh7.n);
%passes_LOOCV8_AUC = sum (passes_LOOCV8) ./ (10001 * stats_mu8.n);
%figure(1); plot (thresholds, passes_LOOCV8,'LineWidth',5);
figure; plot (thresholds, passes_LOOCV7);
% figure; plot (thresholds, passes_LOOCV25,'LineWidth',5);
figure; plot (thresholds, passes_LOOCV3,'LineWidth',5);
% figure; plot (thresholds, passes_LOOCV_Low7);
% figure; plot (thresholds, passes_LOOCV_High7);
figure; plot (thresholds, passes_all);
%figure(3); hist (mu_eff8);
%figure(4); hist (mu_eff7);

temp_paths = Study_paths3;
toss_index_LOOCV3 = find(dice_LOOCV3<0.7);
temp_paths(toss_index_LOOCV3,:) = [];
paths3_pass = temp_paths;
temp_dice = dice_LOOCV3;
temp_dice(toss_index_LOOCV3,:) = [];
dice3_pass = temp_dice;

temp_paths = Study_paths3;
toss_index_LOOCV3 = find(dice_LOOCV3<0.7);
temp_paths(toss_index_LOOCV3,:) = [];
paths3_pass = temp_paths;
temp_dice = dice_LOOCV3;
temp_dice(toss_index_LOOCV3,:) = [];
dice3_pass = temp_dice;

