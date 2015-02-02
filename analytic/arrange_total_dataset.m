function [total,total_all] = arrange_total_dataset ( total );
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;



setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

cd ../../../MATLAB/Tests/direct_search
load total

cd ( path22 );



ix=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1);
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1);
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
total(ix,:) = [];

total_all = total;
if mu_eff_tag(1) == 1    % Eliminate the higher values
    [~, mx] = min( abs( total{1,2}(:,1) - mu_eff_tag(2))); 
    for ii=1:size(total,1)
        aa_iter = total{ii,2};
        aa_iter( mx:end,: ) = []; 
        total{ii,2} = aa_iter;
        aa_iter = total{ii,3};
        aa_iter( mx:end, : ) = [];
        total{ii,3} = aa_iter;
        [total{ii,4}(1), total{ii,4}(3) ] = min(total{ii,2}(:,2));
        total{ii,4}(2) = total{ii,2}( total{ii,4}(3),1);
        [total{ii,5}(1), total{ii,5}(3) ] = max(total{ii,3}(:,7));
        total{ii,5}(2) = total{ii,2}( total{ii,5}(3),1);
    end
    
elseif mu_eff_tag(1) ==2   % Eliminate the lower values
    [~, mx] = min( abs( total{1,2}(:,1) - mu_eff_tag(2))); 
    for ii=1:size(total,1)
        aa_iter = total{ii,2};
        aa_iter( 1:mx, : ) = []; 
        total{ii,2} = aa_iter;
        aa_iter = total{ii,3};
        aa_iter( 1:mx, : ) = [];
        total{ii,3} = aa_iter;
        [total{ii,4}(1), total{ii,4}(3) ] = min(total{ii,2}(:,2));
        total{ii,4}(2) = total{ii,2}( total{ii,4}(3),1);
        [total{ii,5}(1), total{ii,5}(3) ] = max(total{ii,3}(:,7));
        total{ii,5}(2) = total{ii,2}( total{ii,5}(3),1);
    end
end

   
opt.paths = total(:,1);
clear ix mx

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


clear ii jj