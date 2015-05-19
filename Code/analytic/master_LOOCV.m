function [opt, LOOCV, fig_labels] = master_LOOCV ( opttype, data_filename, dice_thresholds, mu_thresholds);
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);

opt.paths = cell (num_studies,2);
for ii = 1:num_studies
    
    opt.paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    opt.paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end

clear ii
indexC = strfind(opt.paths,'Study0035');
toss_index_phantom = find(not(cellfun('isempty',indexC)));
opt.paths(toss_index_phantom,:)=[];
datasummary(toss_index_phantom,:)=[];
num_studies = size(datasummary,1);

indexC = strfind(opt.paths,'0457');
toss_index_phantom = find(not(cellfun('isempty',indexC)));
opt.paths(toss_index_phantom-num_studies,:)=[];
datasummary(toss_index_phantom-num_studies,:)=[];
num_studies = size(datasummary,1);

indexC = strfind(opt.paths,'0476');
toss_index_phantom = find(not(cellfun('isempty',indexC)));
opt.paths(toss_index_phantom-num_studies,:)=[];
datasummary(toss_index_phantom-num_studies,:)=[];
num_studies = size(datasummary,1);

indexC = strfind(opt.paths,'0436');
toss_index_phantom = find(not(cellfun('isempty',indexC)));
opt.paths(toss_index_phantom-num_studies,:)=[];
datasummary(toss_index_phantom-num_studies,:)=[];
num_studies = size(datasummary,1);

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
LOOCV.paths.paths = total_toss;
LOOCV.paths.table   = cell(length_mu_groups, length_dice_thresholds);
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
        
        temp_paths = opt.paths;
        temp_paths(total_toss{ii,jj},:) = [];
        LOOCV.paths.paths{ii,jj} = temp_paths;
        LOOCV.toss_index{ii,jj} = total_toss{ii,jj};

            

        
        if length(total_toss{ii,jj}) < length(opt.paths) -1
            if jj >1
                if length(total_toss{ii,jj}) == length(total_toss{ii,jj-1})
                    
                    opt.mu_eff.values{ii,jj} = opt.mu_eff.values{ii,jj-1};
                    opt.mu_eff.stats{ii,jj}  = opt.mu_eff.stats{ii,jj-1};
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
                    LOOCV.paths.paths{ii,jj} = LOOCV.paths.paths{ii,jj-1};
                    
                    
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
                    opt.mu_eff.stats{ii,jj} = Descriptive_statistics(opt.mu_eff.values{ii,jj});
                    
                    disp(opt.labels{ii,jj});
                    disp(strcat( num2str(kk), [' of '], num2str(length_mu_groups .* length_dice_thresholds), [' groups']));
                    kk = kk+1;
                    [ LOOCV.dice.hh{ii,jj}, LOOCV.dice.values{ii,jj}] = LOOCV_t_test_DiceTemp( LOOCV.paths.paths{ii,jj}, opt.mu_eff.values{ii,jj}, alpha{ii,jj}, best_iter{ii,jj}, opttype);
                    LOOCV.run1(ii,jj) = 2;
                    LOOCV.run2(ii,jj) = 1;
                    LOOCV.dice.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.values{ii,jj});
                    toc
                    
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
                opt.mu_eff.stats{ii,jj} = Descriptive_statistics(opt.mu_eff.values{ii,jj});
                
                disp(opt.labels{ii,jj});
                disp(strcat( num2str(kk), [' of '], num2str(length_mu_groups .* length_dice_thresholds), [' groups']));
                kk = kk+1;
                [ LOOCV.dice.hh{ii,jj}, LOOCV.dice.values{ii,jj}] = LOOCV_t_test_DiceTemp( LOOCV.paths.paths{ii,jj}, opt.mu_eff.values{ii,jj}, alpha{ii,jj}, best_iter{ii,jj}, opttype);
                LOOCV.run1(ii,jj) = 2;
                LOOCV.run2(ii,jj) = 1;
                LOOCV.dice.stats{ii,jj} = Descriptive_statistics_LOOCV( LOOCV.dice.values{ii,jj});
                toc
                
            end
            
            LOOCV.paths.mix{ii,jj}{1,1} = 'Paths';
            LOOCV.paths.mix{ii,jj}{1,2} = 'Dice';
            LOOCV.paths.mix{ii,jj}{1,3} = 'Mu_eff';
            for kk = 2:length(LOOCV.paths.paths{ii,jj})+1
                
                LOOCV.paths.mix{ii,jj}{kk,1} = strcat( LOOCV.paths.paths{ii,jj}{kk-1,1},'/',LOOCV.paths.paths{ii,jj}{kk-1,2});
                LOOCV.paths.mix{ii,jj}{kk,2} = LOOCV.dice.values{ii,jj}(kk-1);
                LOOCV.paths.mix{ii,jj}{kk,3} = opt.mu_eff.values{ii,jj}(kk-1);
            end
            
        else
            kk = kk+1;
            if length(total_toss{ii,jj}) == length(opt.paths) -1
                
                LOOCV.run1(ii,jj) = 1;
                LOOCV.run2(ii,jj) = 0;
            else
                LOOCV.run1(ii,jj) = 0;
                LOOCV.run2(ii,jj) = 0;
            end
        end
        
        
    end
end
toc
end
