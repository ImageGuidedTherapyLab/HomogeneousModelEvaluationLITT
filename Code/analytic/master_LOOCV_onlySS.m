function [opt, LOOCV, fig_labels] = master_LOOCV_onlySS ( total, dice_thresholds, mu_thresholds, naive_mu, opt_tag,summary, choice);

opt.paths = total(:,1);


% variable initialization
if opt_tag ==1
    
    if isempty(dice_thresholds) ==1
        
        opt.labels = 'opt';
        opt.dice.stats = cell(1,1);
        length_dice_thresholds = 0;
    else
        
        % Setup 'total' indexing
        opt_data = cell2mat(total(:,5));
        
        length_dice_thresholds = length(dice_thresholds);
        toss_dice = cell(length_dice_thresholds,1);
        for jj = 1:length_dice_thresholds
            
            toss_dice{jj}= find(dice_thresholds(jj) > opt_data(:,1));
            
        end
        
    end
    clear jj
    
    if isempty(mu_thresholds) ==1
        
        opt.mu_eff.values = cell(1,1);
        length_mu_groups = 1;
        
    else
        
        length_mu_groups = length(mu_thresholds)+1;
        toss_mu = cell(length_mu_groups,1);
        
        
        toss_mu{1} = find( opt_data(:,2) > mu_thresholds(1));
        toss_mu{end} = find(mu_thresholds(end) >= opt_data(:,2));
        
        for ii=2:(length_mu_groups-1)
            
            toss_mu{ii} = find(  (mu_thresholds(ii-1) >= opt_data(:,2)) + (opt_data(:,2) > mu_thresholds(ii))  );
            
        end
    end
    
elseif opt_tag ==2
    if isempty(dice_thresholds) ==1
        
        opt.labels = 'opt';
        opt.dice.stats = cell(1,1);
        length_dice_thresholds = 0;
    else
        
        % Setup 'total' indexing
        opt_data = cell2mat(total(:,4));
        
        for ii = 1:size(opt_data,1)
            dice_data = total{ii,3};
            opt_data(ii,1) = dice_data( opt_data(ii,3), 7);
            
        end
        length_dice_thresholds = length(dice_thresholds);
        toss_dice = cell(length_dice_thresholds,1);
        for jj = 1:length_dice_thresholds
            
            toss_dice{jj}= find(dice_thresholds(jj) > opt_data(:,1));
            
        end
        
    end
    clear jj
    
    if isempty(mu_thresholds) ==1
        
        opt.mu_eff.values = cell(1,1);
        length_mu_groups = 1;
        
    else
        
        length_mu_groups = length(mu_thresholds)+1;
        toss_mu = cell(length_mu_groups,1);
        
        
        toss_mu{1} = find( opt_data(:,2) > mu_thresholds(1));
        toss_mu{end} = find(mu_thresholds(end) >= opt_data(:,2));
        
        for ii=2:(length_mu_groups-1)
            
            toss_mu{ii} = find(  (mu_thresholds(ii-1) >= opt_data(:,2)) + (opt_data(:,2) > mu_thresholds(ii))  );
            
        end
    end
    
end


% length_mu_groups = length(mu_thresholds);
% length_dice_thresholds = length(dice_thresholds);

% Initialize Variables.
total_toss = cell(length_mu_groups, length_dice_thresholds);
opt.mu_eff.values = total_toss;
opt.mu_eff.stats  = total_toss;
opt.mu_eff.run    = total_toss;
opt.mu_eff.run_stat=total_toss;
opt.mu_eff.all.values= opt_data(:,2);
opt.mu_eff.all.stats= Descriptive_statistics(opt.mu_eff.all.values);
opt.dice.values   = total_toss;
opt.dice.stats    = total_toss;

if opt_tag ==1
    opt.dice.all.values = opt_data(:,1);
elseif opt_tag==2
    
    opt.dice.all.values = zeros( size(total_toss,1));
    for ii = 1:size(total,1)
        dice_data = total{ii,3};
        opt.dice.all.values(ii) = dice_data( opt.mu_eff.all.values(ii), 7);
        
    end
    
    
end
opt.dice.all.stats= Descriptive_statistics_LOOCV(opt.dice.all.values);
opt.labels        = total_toss;

opt.dice.all.naive.val = zeros( size(opt_data,1),1);
for ii = 1: size(opt_data,1)
    
    [~, nx] = min( abs( total{1,2}(:,1) - naive_mu)); % naive min
    opt.dice.all.naive.val(ii) = total{ii,3}(nx,7);
end
opt.dice.all.naive.stats = Descriptive_statistics_LOOCV(opt.dice.all.naive.val);

best.group.L2_mean.val       = total_toss;
best.group.L2_mean.stats     = total_toss;
best.group.L2_median.val     = total_toss;
best.group.L2_median.stats   = total_toss;
best.group.dice_mean.val     = total_toss;
best.group.dice_mean.stats   = total_toss;
best.group.dice_median.val   = total_toss;
best.group.dice_median.stats = total_toss;

LOOCV.dice.values         = total_toss;
LOOCV.mu_eff.run          = total_toss;
LOOCV.mu_eff.run_stat     = total_toss;
LOOCV.dice.naive.val      = total_toss;
LOOCV.dice.naive.stats    = total_toss;
LOOCV.dice.stats          = total_toss;
LOOCV.dice.hh             = total_toss;
LOOCV.run1                = zeros(size(total_toss,1),size(total_toss,2));
LOOCV.run2                = zeros(size(total_toss,1),size(total_toss,2));
LOOCV.labels              = total_toss;
LOOCV.paths.paths         = total_toss;
LOOCV.paths.table         = cell(length_mu_groups, length_dice_thresholds);
LOOCV.toss_index          = total_toss;
fig_labels.mu_groups      = cell(size(total_toss,1),1);
fig_labels.DSC            = zeros(size(total_toss,2),1);

best_iter = total_toss;
kk=1;
for ii = 1:length_mu_groups
    for jj = 1:length_dice_thresholds
        
        if isempty(mu_thresholds) ==1   % Identify if there is only one mu_eff group
            total_toss{ii,jj} = toss_dice{jj};
            opt.labels{ii,jj} = strcat(['mu_{eff} all    dice ='], num2str(dice_thresholds(jj)) );
            
        else    % There is >= 1 mu_groups
            total_toss{ii,jj} = unique([ toss_mu{ii}; toss_dice{jj}]);
            if ii ==1    % The first mu_group
                opt.labels{ii,jj} = strcat(['mu{eff} <'], num2str(mu_thresholds(ii)), ['   dice >'], num2str(dice_thresholds(jj)) );
                fig_labels.mu_groups{ii} = strcat(['mu_{eff} <'], num2str(mu_thresholds(ii)));
            elseif ii == length_mu_groups   % The last mu_group
                opt.labels{ii,jj} = strcat(['mu_{eff} >'], num2str(mu_thresholds(ii-1)), ['   dice >'], num2str(dice_thresholds(jj)) );
                fig_labels.mu_groups{ii} = strcat(['mu_{eff} >'], num2str(mu_thresholds(ii-1)));
            else   % Any mu_group that is not first or last (the in-between ones)
                opt.labels{ii,jj} = strcat(['mu_{eff} ='], num2str(mu_thresholds(ii-1)), [' to '], num2str(mu_thresholds(ii)), ['   dice >'], num2str(dice_thresholds(jj)) );
                fig_labels.mu_groups{ii} = strcat(['mu_{eff} ='], num2str(mu_thresholds(ii-1)), [' to '], num2str(mu_thresholds(ii)));
            end
            
        end
        
        fig_labels.DSC(jj) = dice_thresholds(jj);
        
        LOOCV.labels{ii,jj} = opt.labels{ii,jj};
        
        temp_paths = opt.paths;
        temp_paths(total_toss{ii,jj},:) = [];
        LOOCV.paths.paths{ii,jj} = temp_paths;
        LOOCV.toss_index{ii,jj} = total_toss{ii,jj};
        
        if length(total_toss{ii,jj}) < length(opt.paths) -1
            if jj >1
                if length(total_toss{ii,jj}) == length(total_toss{ii,jj-1})
                    
                    opt.mu_eff.values{ii,jj} = opt.mu_eff.values{ii,jj-1};
                    opt.mu_eff.stats{ii,jj}  = opt.mu_eff.stats{ii,jj-1};
                    best_iter{ii,jj} = best_iter{ii,jj-1};
                    opt.dice.values{ii,jj}=opt.dice.values{ii,jj-1};
                    opt.dice.stats{ii,jj} = opt.dice.stats{ii,jj-1};
                    kk=kk+1;
                    LOOCV.dice.hh{ii,jj}=LOOCV.dice.hh{ii,jj-1};
                    LOOCV.run1(ii,jj) = LOOCV.run1(ii,jj-1)+1;
                    LOOCV.run2(ii,jj) = LOOCV.run2(ii,jj-1);
                    LOOCV.dice.stats{ii,jj} = LOOCV.dice.stats{ii,jj-1};
                    LOOCV.dice.values{ii,jj} = LOOCV.dice.values{ii,jj-1};
                    LOOCV.paths.paths{ii,jj} = LOOCV.paths.paths{ii,jj-1};
                    LOOCV.dice.naive.val{ii,jj} = LOOCV.dice.naive.val{ii,jj-1};
                    LOOCV.dice.naive.stats{ii,jj} = LOOCV.dice.naive.stats{ii,jj-1};
                    best.group.L2_mean.val{ii,jj} = best.group.L2_mean.val{ii,jj-1};
                    best.group.L2_mean.stats{ii,jj}     = best.group.L2_mean.stats{ii,jj-1} ;
                    best.group.L2_median.val{ii,jj}     = best.group.L2_median.val{ii,jj-1} ;
                    best.group.L2_median.stats{ii,jj}   = best.group.L2_median.stats{ii,jj-1};
                    best.group.dice_mean.val{ii,jj}     = best.group.dice_mean.val{ii,jj-1} ;
                    best.group.dice_mean.stats{ii,jj}   = best.group.dice_mean.stats{ii,jj-1};
                    best.group.dice_median.val{ii,jj}   = best.group.dice_median.val{ii,jj-1};
                    best.group.dice_median.stats{ii,jj} =  best.group.dice_median.stats{ii,jj-1};
                    
                    
                else
                    total_iter = total;
                    total_iter( total_toss{ii,jj} , :) = [];
                    
                    num_runs = size(total_iter,1);
                    num_mu   = size(total_iter{1,2},1);
                    L2 = zeros( num_mu, num_runs);
                    dice = L2;
                    for mm = 1:num_runs
                        
                        L2(:,mm) = total_iter{mm,2}(:,2);
                        dice(:,mm) = total_iter{mm,3}(:,7);
                        
                    end
                    clear mm
                    
                    L2_mean = mean( L2,2);
                    L2_median = median( L2,2);
                    dice_mean = mean( dice,2);
                    dice_median = median( dice,2);
                    
                    [~, L2_mean_ix]= min( L2_mean);
                    [~, L2_median_ix] = min( L2_median );
                    [~, dice_mean_ix] = max( dice_mean );
                    [~, dice_median_ix] = max( dice_median );
                    
                    opt.mu_eff.values{ii,jj} = opt_data(:,2);
                    opt.mu_eff.values{ii,jj} (total_toss{ii,jj})=[];
                    
                    best_iter{ii,jj} = opt_data(:,3);
                    best_iter{ii,jj} (total_toss{ii,jj})=[];
                    
                    opt.dice.values{ii,jj} = opt_data(:,1);
                    opt.dice.values{ii,jj}(total_toss{ii,jj}) = [];
                    
                    opt.dice.stats{ii,jj} = Descriptive_statistics_LOOCV(opt.dice.values{ii,jj});
                    opt.mu_eff.stats{ii,jj} = Descriptive_statistics(opt.mu_eff.values{ii,jj});
                    
                    disp(opt.labels{ii,jj});
                    disp(strcat( num2str(kk), [' of '], num2str(length_mu_groups .* length_dice_thresholds), [' groups']));
                    kk = kk+1;
                    
                    mu_iter = opt.mu_eff.values{ii,jj};
                    length_iter = length(opt.mu_eff.values{ii,jj});
                    LOOCV.dice.values{ii,jj} = zeros(length_iter,1);
                    LOOCV.dice.naive.val{ii,jj} = zeros(length_iter,1);
                    LOOCV.mu_eff.run{ii,jj}  = zeros(length_iter,1);
                    best.group.L2_mean.val{ii,jj}=zeros(length_iter,1);
                    best.group.L2_median.val{ii,jj}  = zeros(length_iter,1);
                    best.group.dice_mean.val{ii,jj}  = zeros(length_iter,1);
                    best.group.dice_median.val{ii,jj}= zeros(length_iter,1);
                    
                    for ll = 1:length_iter
                        
                        mu = mu_iter;
                        mu(ll) = [];
                        mu = round(mean(mu));
                        LOOCV.mu_eff.run{ii,jj}(ll)=mu;
                        [~, ix] = min( abs( total_iter{1,2}(:,1) - mu));
                        [~, nx] = min( abs( total_iter{1,2}(:,1) - naive_mu)); % naive min
                        LOOCV.dice.values{ii,jj}(ll) = total_iter{ll,3}(ix,7);
                        LOOCV.dice.naive.val{ii,jj}(ll) = total_iter{ll,3}(nx,7);
                        
                        [~, ix] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( L2_mean_ix,1)));
                        [~, jx] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( L2_median_ix,1))); % naive min
                        [~, kx] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( dice_mean_ix,1)));
                        [~, lx] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( dice_median_ix,1))); % naive min
                        best.group.L2_mean.val{ii,jj}(ll)        = total_iter{ll,3}(ix,7);
                        best.group.L2_median.val{ii,jj}(ll)      = total_iter{ll,3}(jx,7);
                        best.group.dice_mean.val{ii,jj}(ll)      = total_iter{ll,3}(kx,7);
                        best.group.dice_median.val{ii,jj}(ll)    = total_iter{ll,3}(lx,7);
                    end
                    [ LOOCV.dice.hh{ii,jj}.H,  LOOCV.dice.hh{ii,jj}.ptest,  LOOCV.dice.hh{ii,jj}.ci,  LOOCV.dice.hh{ii,jj}.stats] = ttest( LOOCV.dice.values{ii,jj}, 0.7, 0.05, 'right');
                    LOOCV.run1(ii,jj) = 2;
                    LOOCV.run2(ii,jj) = 1;
                    LOOCV.dice.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.values{ii,jj});
                    LOOCV.dice.naive.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.naive.val{ii,jj});
                    LOOCV.mu_eff.run_stat{ii,jj} = Descriptive_statistics( LOOCV.mu_eff.run{ii,jj});
                    best.group.L2_mean.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.L2_mean.val{ii,jj});
                    best.group.L2_median.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.L2_median.val{ii,jj});
                    best.group.dice_mean.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.dice_mean.val{ii,jj});
                    best.group.dice_median.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.dice_median.val{ii,jj});
                    
                end
            else
                
                total_iter = total;
                total_iter( total_toss{ii,jj} , :) = [];
                
                num_runs = size(total_iter,1);
                num_mu   = size(total_iter{1,2},1);
                L2 = zeros( num_mu, num_runs);
                dice = L2;
                for mm = 1:num_runs
                    
                    L2(:,mm) = total_iter{mm,2}(:,2);
                    dice(:,mm) = total_iter{mm,3}(:,7);
                    
                end
                clear mm
                
                L2_mean = mean( L2,2);
                L2_median = median( L2,2);
                dice_mean = mean( dice,2);
                dice_median = median( dice,2);
                
                [~, L2_mean_ix]= min( L2_mean);
                [~, L2_median_ix] = min( L2_median );
                [~, dice_mean_ix] = max( dice_mean );
                [~, dice_median_ix] = max( dice_median );
                
                opt.mu_eff.values{ii,jj} = opt_data(:,2);
                opt.mu_eff.values{ii,jj} (total_toss{ii,jj})=[];
                
                best_iter{ii,jj} = opt_data(:,3);
                best_iter{ii,jj} (total_toss{ii,jj})=[];
                
                opt.dice.values{ii,jj} = opt_data(:,1);
                opt.dice.values{ii,jj}(total_toss{ii,jj}) = [];
                
                opt.dice.stats{ii,jj} = Descriptive_statistics_LOOCV(opt.dice.values{ii,jj});
                opt.mu_eff.stats{ii,jj} = Descriptive_statistics(opt.mu_eff.values{ii,jj});
                
                disp(opt.labels{ii,jj});
                disp(strcat( num2str(kk), [' of '], num2str(length_mu_groups .* length_dice_thresholds), [' groups']));
                kk = kk+1;
                
                mu_iter = opt.mu_eff.values{ii,jj};
                length_iter = length(opt.mu_eff.values{ii,jj});
                LOOCV.dice.values{ii,jj} = zeros(length_iter,1);
                LOOCV.dice.naive.val{ii,jj} = zeros(length_iter,1);
                LOOCV.mu_eff.run{ii,jj}  = zeros(length_iter,1);
                best.group.L2_mean.val{ii,jj}=zeros(length_iter,1);
                best.group.L2_median.val{ii,jj}  = zeros(length_iter,1);
                best.group.dice_mean.val{ii,jj}  = zeros(length_iter,1);
                best.group.dice_median.val{ii,jj}= zeros(length_iter,1);
                
                for ll = 1:length_iter
                    
                    mu = mu_iter;
                    mu(ll) = [];
                    mu = round(mean(mu));
                    LOOCV.mu_eff.run{ii,jj}(ll)=mu;
                    [~, ix] = min( abs( total_iter{1,2}(:,1) - mu));
                    [~, nx] = min( abs( total_iter{1,2}(:,1) - naive_mu)); % naive min
                    LOOCV.dice.values{ii,jj}(ll) = total_iter{ll,3}(ix,7);
                    LOOCV.dice.naive.val{ii,jj}(ll) = total_iter{ll,3}(nx,7);
                    
                    [~, ix] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( L2_mean_ix,1)));
                    [~, jx] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( L2_median_ix,1))); % naive min
                    [~, kx] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( dice_mean_ix,1)));
                    [~, lx] = min( abs( total_iter{1,2}(:,1) - total_iter{1,2}( dice_median_ix,1))); % naive min
                    best.group.L2_mean.val{ii,jj}(ll)        = total_iter{ll,3}(ix,7);
                    best.group.L2_median.val{ii,jj}(ll)      = total_iter{ll,3}(jx,7);
                    best.group.dice_mean.val{ii,jj}(ll)      = total_iter{ll,3}(kx,7);
                    best.group.dice_median.val{ii,jj}(ll)    = total_iter{ll,3}(lx,7);
                end
                [ LOOCV.dice.hh{ii,jj}.H,  LOOCV.dice.hh{ii,jj}.ptest,  LOOCV.dice.hh{ii,jj}.ci,  LOOCV.dice.hh{ii,jj}.stats] = ttest( LOOCV.dice.values{ii,jj}, 0.7, 0.05, 'right');
                LOOCV.run1(ii,jj) = 2;
                LOOCV.run2(ii,jj) = 1;
                LOOCV.dice.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.values{ii,jj});
                LOOCV.dice.naive.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.naive.val{ii,jj});
                LOOCV.mu_eff.run_stat{ii,jj} = Descriptive_statistics( LOOCV.mu_eff.run{ii,jj});
                best.group.L2_mean.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.L2_mean.val{ii,jj});
                best.group.L2_median.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.L2_median.val{ii,jj});
                best.group.dice_mean.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.dice_mean.val{ii,jj});
                best.group.dice_median.stats{ii,jj} = Descriptive_statistics_LOOCV( best.group.dice_median.val{ii,jj});
                
            end
            LOOCV.paths.mix{ii,jj} = cell( length(LOOCV.paths.paths{ii,jj}),3);
            LOOCV.paths.mix{ii,jj}{1,1} = 'Paths';
            LOOCV.paths.mix{ii,jj}{1,2} = 'Dice';
            LOOCV.paths.mix{ii,jj}{1,3} = 'Mu_eff';
            for kk = 1:length(LOOCV.paths.paths{ii,jj})
                %  LOOCV.paths.mix{ii,jj}{kk,1} = strcat( LOOCV.paths.paths{ii,jj}{kk-1,1},'/',LOOCV.paths.paths{ii,jj}{kk-1,2});
                LOOCV.paths.mix{ii,jj}{kk+1,1} = LOOCV.paths.paths{ii,jj}{kk,1};
                LOOCV.paths.mix{ii,jj}{kk+1,2} = LOOCV.dice.values{ii,jj}(kk);
                LOOCV.paths.mix{ii,jj}{kk+1,3} = opt.mu_eff.values{ii,jj}(kk);
            end
            
        else
            kk = kk+1;
            if length(total_toss{ii,jj}) == length(opt.paths) -1
                
                LOOCV.run1(ii,jj) = 1;  % Indicates only one dataset qualifies for this group
                LOOCV.run2(ii,jj) = 0;  % Run canceled
            else
                LOOCV.run1(ii,jj) = 0; % Indicates no dataset qualifies for this group
                LOOCV.run2(ii,jj) = 0; % Run canceled
            end
        end
        
        
    end
end

end
