function [opt, LOOCV, fig_labels] = master_LOOCV ( opttype, data_filename, dice_thresholds, mu_thresholds, Matlab_flag);

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
    
    opt.labels = 'opt';
    opt.dice.stats = cell(1,1);
    length_dice_thresholds = 0;
else
    
    length_dice_thresholds = length(dice_thresholds);
    toss_dice = cell(length_dice_thresholds,1);
    for jj = 1:length_dice_thresholds
        
        toss_dice{jj}= find(dice_thresholds(jj) > datasummary(:,7));
        
    end
    
end
clear jj
if isempty(mu_thresholds) ==1
    
    opt.mu_eff.values = cell(1,1);
    length_mu_groups = 1;
    
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

% Initialize Variables.
total_toss = cell(length_mu_groups, length_dice_thresholds);
opt.mu_eff.values = total_toss;
opt.mu_eff.stats  = total_toss;
opt.mu_eff.all.values= datasummary(:,4);
opt.mu_eff.all.stats= Descriptive_statistics(opt.mu_eff.all.values);
opt.dice.values   = total_toss;
opt.dice.stats    = total_toss;
opt.dice.all.values = datasummary(:,7);
opt.dice.all.stats= Descriptive_statistics(opt.dice.all.values);
opt.labels        = total_toss;

LOOCV.dice.values = total_toss;
LOOCV.dice.stats  = total_toss;
LOOCV.dice.hh     = total_toss;
LOOCV.run1        = zeros(size(total_toss,1),size(total_toss,2));
LOOCV.run2        = zeros(size(total_toss,1),size(total_toss,2));
LOOCV.labels      = total_toss;
LOOCV.paths       = total_toss;
LOOCV.toss_index  = total_toss;
fig_labels.mu_groups = cell(size(total_toss,1),1);
fig_labels.DSC    = zeros(size(total_toss,2),1);

alpha = total_toss;
best_iter = total_toss;
kk=1;
for ii = 1:length_mu_groups
    for jj = 1:length_dice_thresholds
        
        if isempty(mu_thresholds) ==1   % Identify if there is only one mu_eff group
            total_toss{ii,jj} = toss_dice{jj};
            opt.labels{ii,jj} = strcat(['mu_eff all    dice ='], num2str(dice_thresholds(jj)) );
            
        else    % There is >= 1 mu_groups
            total_toss{ii,jj} = unique([ toss_mu{ii}; toss_dice{jj}]);
            if ii ==1    % The first mu_group
                opt.labels{ii,jj} = strcat(['mu_eff <'], num2str(mu_thresholds(ii)), ['   dice >'], num2str(dice_thresholds(jj)) );
                fig_labels.mu_groups{ii} = strcat(['mu_eff <'], num2str(mu_thresholds(ii)));
            elseif ii == length_mu_groups   % The last mu_group
                opt.labels{ii,jj} = strcat(['mu_eff >'], num2str(mu_thresholds(ii-1)), ['   dice >'], num2str(dice_thresholds(jj)) );
                fig_labels.mu_groups{ii} = strcat(['mu_eff >'], num2str(mu_thresholds(ii-1)));
            else   % Any mu_group that is not first or last (the in-between ones)
                opt.labels{ii,jj} = strcat(['mu_eff ='], num2str(mu_thresholds(ii-1)), [' to '], num2str(mu_thresholds(ii)), ['   dice >'], num2str(dice_thresholds(jj)) );
                fig_labels.mu_groups{ii} = strcat(['mu_eff ='], num2str(mu_thresholds(ii-1)), [' to '], num2str(mu_thresholds(ii)));
            end
            
        end
        
        fig_labels.DSC(jj) = dice_thresholds(jj); 
        
        LOOCV.labels{ii,jj} = opt.labels{ii,jj};
        
        temp_paths = Study_paths;
        temp_paths(total_toss{ii,jj},:) = [];
        LOOCV.paths{ii,jj} = temp_paths;
        LOOCV.toss_index{ii,jj} = total_toss{ii,jj};
        
        if length(total_toss{ii,jj}) < length(Study_paths) -1
            if jj >1
                if length(total_toss{ii,jj}) == length(total_toss{ii,jj-1})
                    
                    opt.mu_eff.values{ii,jj} = opt.mu_eff.values{ii,jj-1};
                    alpha{ii,jj} = alpha{ii,jj-1};
                    best_iter{ii,jj} = best_iter{ii,jj-1};
                    opt.dice.values{ii,jj}=opt.dice.values{ii,jj-1};
                    opt.dice.stats{ii,jj} = opt.dice.stats{ii,jj-1};
                    kk=kk+1;
                    LOOCV.dice.hh{ii,jj}=LOOCV.dice.hh{ii,jj-1};
                    LOOCV.run1(ii,jj) = LOOCV.run1(ii,jj-1)+1;
                    LOOCV.run2(ii,jj) = LOOCV.run2(ii,jj-1);
                    LOOCV.dice.stats{ii,jj} = LOOCV.dice.stats{ii,jj-1};
                    LOOCV.dice.values{ii,jj} = LOOCV.dice.values{ii,jj-1};
                    
                else
                    opt.mu_eff.values{ii,jj} = datasummary(:,4);
                    opt.mu_eff.values{ii,jj} (total_toss{ii,jj})=[];
                    
                    alpha{ii,jj} = datasummary(:,5);
                    alpha{ii,jj} (total_toss{ii,jj}) = [];
                    
                    best_iter{ii,jj} = datasummary(:,3);
                    best_iter{ii,jj} (total_toss{ii,jj})=[];
                    
                    opt.dice.values{ii,jj} = datasummary(:,7);
                    opt.dice.values{ii,jj}(total_toss{ii,jj}) = [];
                    
                    opt.dice.stats{ii,jj} = Descriptive_statistics(opt.dice.values{ii,jj});
                    opt.dice.stats{ii,jj} = Descriptive_statistics(opt.mu_eff.values{ii,jj});
                    
                    disp(opt.labels{ii,jj});
                    disp(strcat( num2str(kk), [' of '], num2str(length_mu_groups .* length_dice_thresholds), [' groups']));
                    kk = kk+1;
                    [ LOOCV.dice.hh{ii,jj}, LOOCV.dice.values{ii,jj}] = LOOCV_t_test_DiceTemp( LOOCV.paths{ii,jj}, opt.mu_eff.values{ii,jj}, alpha{ii,jj}, best_iter{ii,jj}, opttype, Matlab_flag);
                    LOOCV.run1(ii,jj) = 2;
                    LOOCV.run2(ii,jj) = 1;
                    LOOCV.dice.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.values{ii,jj});
                    
                end
            else
                
                opt.mu_eff.values{ii,jj} = datasummary(:,4);
                opt.mu_eff.values{ii,jj} (total_toss{ii,jj})=[];
                
                alpha{ii,jj} = datasummary(:,5);
                alpha{ii,jj} (total_toss{ii,jj}) = [];
                
                best_iter{ii,jj} = datasummary(:,3);
                best_iter{ii,jj} (total_toss{ii,jj})=[];
                
                opt.dice.values{ii,jj} = datasummary(:,7);
                opt.dice.values{ii,jj}(total_toss{ii,jj}) = [];
                
                opt.dice.stats{ii,jj} = Descriptive_statistics(opt.dice.values{ii,jj});
                opt.dice.stats{ii,jj} = Descriptive_statistics(opt.mu_eff.values{ii,jj});
                
                disp(opt.labels{ii,jj});
                disp(strcat( num2str(kk), [' of '], num2str(length_mu_groups .* length_dice_thresholds), [' groups']));
                kk = kk+1;
                [ LOOCV.dice.hh{ii,jj}, LOOCV.dice.values{ii,jj}] = LOOCV_t_test_DiceTemp( LOOCV.paths{ii,jj}, opt.mu_eff.values{ii,jj}, alpha{ii,jj}, best_iter{ii,jj}, opttype, Matlab_flag);
                LOOCV.run1(ii,jj) = 2;
                LOOCV.run2(ii,jj) = 1;
                LOOCV.dice.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.values{ii,jj});
            end
            
        else
            kk = kk+1;
            if length(total_toss{ii,jj}) == length(Study_paths) -1
                
                LOOCV.run1(ii,jj) = 1;
                LOOCV.run2(ii,jj) = 0;
            else
                LOOCV.run1(ii,jj) = 0;
                LOOCV.run2(ii,jj) = 0;
            end
        end
        
        
    end
end

end

%
%
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
% %figure(3); hist (mu_eff8);
% %figure(4); hist (mu_eff7);
%
% temp_paths = Study_paths3;
% toss_index_LOOCV3 = find(dice_LOOCV3<0.7);
% temp_paths(toss_index_LOOCV3,:) = [];
% paths3_pass = temp_paths;
% temp_dice = dice_LOOCV3;
% temp_dice(toss_index_LOOCV3,:) = [];
% dice3_pass = temp_dice;
%
% temp_paths = Study_paths3;
% toss_index_LOOCV3 = find(dice_LOOCV3<0.7);
% temp_paths(toss_index_LOOCV3,:) = [];
% paths3_pass = temp_paths;
% temp_dice = dice_LOOCV3;
% temp_dice(toss_index_LOOCV3,:) = [];
% dice3_pass = temp_dice;

