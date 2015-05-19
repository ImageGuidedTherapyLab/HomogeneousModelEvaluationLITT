function [opt, LOOCV, fig_labels] = master_LOOCV_onlySS_rand ( total, dice_thresholds, naive_var, opt_tag,choice);
tic
opt.paths = total(:,1);

 cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
% data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file
% 
opttype = 'bestfit50';

[~, ~,max_phys_sz] =  display_inputvars ( 1:30);


% datasummary = dlmread(data_filename,',',1,0);
% datasummary(any(isnan(datasummary), 2), 7) = 1;
% num_studies = size(datasummary,1);
% for ii = 1:num_studies
%     
%     Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
%     Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
%     
% end
% clear ii
% 
% Study_paths(toss_ix,:)=[];


% variable initialization
if opt_tag ==1
    
    if isempty(dice_thresholds) ==1
        
        opt.labels = 'opt';
        opt.dice.stats = cell(1,1);
        length_dice_thresholds = 0;
    else
        
        % Setup 'total' indexing
        opt_data = cell2mat(total(:,8));
        
        length_dice_thresholds = length(dice_thresholds);
        toss_dice = cell(length_dice_thresholds,1);
        for jj = 1:length_dice_thresholds
            
            toss_dice{jj}= find(dice_thresholds(jj) > opt_data(:,1));
            
        end
        
    end
    clear jj
    
elseif opt_tag ==2  % optimize on L2
    if isempty(dice_thresholds) ==1
        
        opt.labels = 'opt';
        opt.dice.stats = cell(1,1);
        length_dice_thresholds = 0;
    else
        
        % Setup 'total' indexing
        aa = cell2mat( total(:,7) );  % Find the optimal L2 info
        sz_a = size(aa);
        opt_L2 = aa( (sz_a(1)/3 +1): (sz_a(1)/3*2),:);    % Keep the L2 portion
        opt_data = zeros( size(opt_L2,1),2);
        opt_data(:,2) = opt_L2(:,3);
        for ii = 1:size( opt_L2,1)
            
            opt_data(ii,1) = total{ii,3}( opt_data(ii,2),7);  % Finds the DSC corresponding to the L2
            
        end
        
        length_dice_thresholds = length(dice_thresholds);
        toss_dice = cell(length_dice_thresholds,1);
        for jj = 1:length_dice_thresholds
            
            toss_dice{jj}= find(dice_thresholds(jj) > opt_data(:,1));
            
        end
        
    end
    clear jj
    
    
end


% length_dice_thresholds = length(dice_thresholds);

% Initialize Variables.
total_toss = cell(length_dice_thresholds, 1);
opt.var.values = total_toss;
opt.var.stats.mu  = total_toss;
opt.var.stats.perf  = total_toss;
opt.var.run    = total_toss;
opt.var.run_stat=total_toss;
opt.var.all.values= opt_data(:,2:3);
opt.var.all.stats.mu= Descriptive_statistics(opt.var.all.values(:,1));
opt.var.all.stats.perf= Descriptive_statistics(opt.var.all.values(:,2));
opt.dice.values   = total_toss;
opt.dice.stats    = total_toss;

if opt_tag ==1  % optimize on 57 C isotherm DSC
    opt.dice.all.values = opt_data(:,1);
elseif opt_tag==2  % optimize on L2
    
    opt.dice.all.values = opt_data(:,1);
    
end
opt.dice.all.stats= Descriptive_statistics_LOOCV(opt.dice.all.values);
opt.labels        = total_toss;

opt.dice.all.naive.val = zeros( size(opt_data,1),1);
for ii = 1: size(opt_data,1)
    
    [~, nx] = min( abs( total{1,2}(:,1) - naive_var)); % naive min
    opt.dice.all.naive.val(ii) = total{ii,3}(nx,7);
end
opt.dice.all.naive.stats = Descriptive_statistics_LOOCV(opt.dice.all.naive.val);

% best.group.L2_mean.val       = total_toss;
% best.group.L2_mean.stats     = total_toss;
% best.group.L2_median.val     = total_toss;
% best.group.L2_median.stats   = total_toss;
% best.group.dice_mean.val     = total_toss;
% best.group.dice_mean.stats   = total_toss;
% best.group.dice_median.val   = total_toss;
% best.group.dice_median.stats = total_toss;

LOOCV.dice.values         = total_toss;
LOOCV.var.run          = total_toss;
LOOCV.var.run_stat     = total_toss;
LOOCV.dice.naive.val      = total_toss;
LOOCV.dice.naive.stats    = total_toss;
LOOCV.dice.stats          = total_toss;
LOOCV.dice.hh             = total_toss;
LOOCV.run1                = zeros(size(total_toss,1),size(total_toss,2));
LOOCV.run2                = zeros(size(total_toss,1),size(total_toss,2));
LOOCV.labels              = total_toss;
LOOCV.paths.paths         = total_toss;
LOOCV.paths.table         = cell( length_dice_thresholds,1);
LOOCV.toss_index          = total_toss;
fig_labels.var_groups      = cell(size(total_toss,1),1);
fig_labels.DSC            = zeros(size(total_toss,2),1);

%best_iter = total_toss;
kk=1;
for ii = 1:length_dice_thresholds
    
    
    
    total_toss{ii} = toss_dice{ii};
    opt.labels{ii} = strcat(['Var all    dice ='], num2str(dice_thresholds(ii)) );
    
    
    
    fig_labels.DSC(ii) = dice_thresholds(ii);
    
    LOOCV.labels{ii} = opt.labels{ii};
    
    temp_paths = opt.paths;
    temp_paths(total_toss{ii},:) = [];
    LOOCV.paths.paths{ii} = temp_paths;
    LOOCV.toss_index{ii} = total_toss{ii};
    
    if length(total_toss{ii}) < length(opt.paths) -1
        
        
        
        total_iter = total;
        total_iter( total_toss{ii} , :) = [];
        
        %                 num_runs = size(total_iter,1);
        %                 num_var   = size(total_iter{1,2},1);
        %                 L2 = zeros( num_var, num_runs);
        %                 dice = L2;
        %                 for mm = 1:num_runs
        %
        %                     L2(:,mm) = total_iter{mm,2}(:,3);
        %                     dice(:,mm) = total_iter{mm,3}(:,7);
        %
        %                 end
        %                 clear mm
        %
        %                 L2_mean = mean( L2,2);
        %                 L2_median = median( L2,2);
        %                 dice_mean = mean( dice,2);
        %                 dice_median = median( dice,2);
        %
        %                 [~, L2_mean_ix]= min( L2_mean);
        %                 [~, L2_median_ix] = min( L2_median );
        %                 [~, dice_mean_ix] = max( dice_mean );
        %                 [~, dice_median_ix] = max( dice_median );
        
        opt.var.values{ii} = opt_data(:,2:3);
        opt.var.values{ii} (total_toss{ii})=[];
        
        %                 best_iter{ii,jj} = opt_data(:,3);
        %                 best_iter{ii,jj} (total_toss{ii,jj})=[];
        
        opt.dice.values{ii} = opt_data(:,1);
        opt.dice.values{ii}(total_toss{ii}) = [];
        
        opt.dice.stats{ii} = Descriptive_statistics_LOOCV(opt.dice.values{ii});
        opt.var.stats.mu{ii} = Descriptive_statistics(opt.var.values{ii}(:,1));
        opt.var.stats.perf{ii} = Descriptive_statistics(opt.var.values{ii}(:,2));
        
        disp(opt.labels{ii});
        disp(strcat( num2str(kk), [' of '], num2str( length_dice_thresholds), [' groups']));
        kk = kk+1;
        
        var_iter = opt.var.values{ii};
        length_iter = length(opt.var.values{ii});
        LOOCV.dice.values{ii} = zeros(length_iter,1);
        LOOCV.dice.naive.val{ii} = zeros(length_iter,2);
        LOOCV.var.run{ii}  = zeros(length_iter,2);
        %                 best.group.L2_mean.val{ii,jj}=zeros(length_iter,1);
        %                 best.group.L2_median.val{ii,jj}  = zeros(length_iter,1);
        %                 best.group.dice_mean.val{ii,jj}  = zeros(length_iter,1);
        %                 best.group.dice_median.val{ii,jj}= zeros(length_iter,1);
        
        for ll = 1:length_iter
            
            var = var_iter;
            var(ll,:) = [];
            if choice ==4
                var = median(var,1);
                LOOCV.var.run{ii}(ll,:)=var;
                
            elseif choice==5
            var = mean(var,1);
            LOOCV.var.run{ii}(ll,:)=var;
            
            path_base = strcat ( 'workdir/',total_iter{ll,1}, '/opt');
            load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );
            
            %[all_opt_fig, no_pwr_fig,sim_dim] = temperature_GPU_LOOCV_rerun ( inputdatavars, 50, max_phys_sz, mu_eff_list, w_perf, k_cond, choice );
            [LOOCV.dice.values{ii}(ll)] = temperature_GPU_LOOCV_rerun ( inputdatavars, 50, max_phys_sz, var(1), var(2), 0.527, 1 );
            end
            
%             [~, ix] = min( abs( total_iter{1,2}(:,1) - var));
%             [~, iix] = min( abs( total_iter{1,2}(:,1) - var));
%             [~, nx] = min( abs( total_iter{1,2}(:,1) - naive_var)); % naive min
%             LOOCV.dice.values{ii}(ll) = total_iter{ll,3}(ix,7);
%             LOOCV.dice.naive.val{ii}(ll) = total_iter{ll,3}(nx,7);
            
        end
        [ LOOCV.dice.hh{ii}.H,  LOOCV.dice.hh{ii}.ptest,  LOOCV.dice.hh{ii}.ci,  LOOCV.dice.hh{ii}.stats] = ttest( LOOCV.dice.values{ii}, 0.7, 0.05, 'right');
        LOOCV.run1(ii) = 2;
        LOOCV.run2(ii) = 1;
        LOOCV.dice.stats{ii} = Descriptive_statistics_LOOCV( LOOCV.dice.values{ii});
        LOOCV.dice.naive.stats{ii} = Descriptive_statistics_LOOCV( LOOCV.dice.naive.val{ii});
        LOOCV.var.run_stat{ii} = Descriptive_statistics( LOOCV.var.run{ii});
        %                 best.group.L2_mean.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.L2_mean.val{ii,jj});
        %                 best.group.L2_median.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.L2_median.val{ii,jj});
        %                 best.group.dice_mean.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.dice_mean.val{ii,jj});
        %                 best.group.dice_median.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.dice_median.val{ii,jj});
        
        
        LOOCV.paths.mix{ii} = cell( length(LOOCV.paths.paths{ii}),3);
        LOOCV.paths.mix{ii}{1,1} = 'Paths';
        LOOCV.paths.mix{ii}{1,2} = 'Dice';
        LOOCV.paths.mix{ii}{1,3} = 'Var';
        for kk = 1:length(LOOCV.paths.paths{ii})
            %  LOOCV.paths.mix{ii,jj}{kk,1} = strcat( LOOCV.paths.paths{ii,jj}{kk-1,1},'/',LOOCV.paths.paths{ii,jj}{kk-1,2});
            LOOCV.paths.mix{ii}{kk+1,1} = LOOCV.paths.paths{ii}{kk,1};
            LOOCV.paths.mix{ii}{kk+1,2} = LOOCV.dice.values{ii}(kk);
            LOOCV.paths.mix{ii}{kk+1,3} = opt.var.values{ii}(kk);
        end
        
    else
        kk = kk+1;
        if length(total_toss{ii}) == length(opt.paths) -1
            
            LOOCV.run1(ii) = 1;  % Indicates only one dataset qualifies for this group
            LOOCV.run2(ii) = 0;  % Run canceled
        else
            LOOCV.run1(ii) = 0; % Indicates no dataset qualifies for this group
            LOOCV.run2(ii) = 0; % Run canceled
        end
    end
    
    
    
end

end
