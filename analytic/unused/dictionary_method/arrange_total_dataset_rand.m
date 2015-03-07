%function [total,best,total_all] = arrange_total_dataset ( data_filename, var_tag, choice );

function [total,total_all,summary,ix] = arrange_total_dataset_rand ( data_filename, var_tag,toss_choice );
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

%cd ../../../MATLAB/Tests/direct_search
cd ../../../MATLAB/Tests/direct_search/libraries

load ('GPU_dict_perf_mu_rand.mat');

cd ( path22 );

% List of excluded datasets
total(1,:) = []; % Drop the labels
ix = 0;
if toss_choice == 1
    ix       =find(~cellfun(@isempty,regexp(total(:,1),'0497'))==1); % Absolutely should be excluded
    %total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0378'))==1); % Strongly suggest exclusion
    %total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1); % Strongly suggest exclusion
    %total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1); % Absolutely should be excluded
    %total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0466'))==1); % Very probably suggest exclusion
    %total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0468'))==1); % Very probably suggest exclusion
    %total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0471'))==1); % Strongly suggest exclusion
    %total(ix,:) = [];
    %
    % ix=find(~cellfun(@isempty,regexp(total(:,1),'0417'))==1); % Very probably suggest exclusion
    % total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0409'))==1); % Absolutely should be excluded
    %total(ix,:) = [];
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0415'))==1); % Absolutely should be excluded
    total(ix,:) = [];
end


% ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
% total(ix,:) = [];

total_all = total;
if var_tag(1,2) == 1
    ix = find( total{1,2}(:,1) > 300);
elseif var_tag(1,2) == 2
    ix = find( total{1,2}(:,1) <1500);
end

if var_tag(1,2) == 1 || var_tag(1,2)==2
    for ii=1:size(total,1)
        
        ix = find( total{ii,2}(:,1) > 300);
        for kk=2:5
            total{ii,kk}(ix,:) = [];
        end
        total{ii,6}(ix,:,:)=[];
        
        %         aa_iter = total{ii,2};
        %         aa_iter( mx:end,: ) = [];
        %         total{ii,2} = aa_iter;
        %         aa_iter = total{ii,3};
        %         aa_iter( mx:end, : ) = [];
        %         total{ii,3} = aa_iter;
        %         aa_iter = total{ii,4};
        %         aa_iter( mx:end, : else
  
        %         total{ii,4} = aa_iter;
        %         aa_iter = total{ii,5};
        %         aa_iter( mx:end, : ) = [];
        %         total{ii,5} = aa_iter;
        %         aa_iter = total{ii,6};
        %         aa_iter(mx:end, :, :)= [];
        %         total{ii,6} = aa_iter;
        
        [total] = optimal_metrics_rand (total,ii);

% Record optimal L information
% total{ii,7} = zeros(3,4);
% [total{ii,7}(1,1) , index] = min (total{ii,2}(:,3));  % record optimal L1
% total{ii,7}(1,2:3) = total{ii,2}(index,1:2);   % record values that produces optimal L1
% total{ii,7}(1,4) = index; % record index that produces optimal L1
% [total{ii,7}(2,1) , index] = min (total{ii,2}(:,4));  % record optimal L2
% total{ii,7}(2,2:3) = total{ii,2}(index,1:2);   % record value that produces optimal L2
% total{ii,7}(2,4) = index; % record index that produces optimal L2
% [total{ii,7}(3,1) , index] = min (total{ii,2}(:,5));  % record optimal L_inf
% total{ii,7}(3,2:3) = total{ii,2}(index,1:2);   % record value that produces optimal L_inf
% total{ii,7}(3,4) = index; % record index that produces optimal L_inf
% 
% % Record optimal 57 C isotherm DSC information
% total{ii,8} = zeros(1,4);
% [total{ii,8}(1) , index] = max (total{ii,3}(:,7));  % record optimal Dice
% total{ii,8}(2:3) = total{ii,2}(index,1:2);   % record value that produces optimal Dice
% total{ii,8}(4) = index; % record index that produces optimal Dice
% 
% % Record optimal 57 C Hausdorff distance information
% total{ii,9} = zeros(1,4);
% [total{ii,9}(1) , index] = min (total{ii,4}(:,7)); % record optimal Hausdorff Distance
% total{ii,9}(2:3) = total{ii,2}(index,1:2); % record value that produces optimal Hausdorff distance
% total{ii,9}(4) = index;
% 
% % Record optimal temperature MI and 57 C isotherm MI
% total{ii,10} = zeros(2,4);
% [total{ii,10}(1,1) , index] = max (total{ii,2}(:,5));  % MI for temperature
% total{ii,10}(1,2:3) = total{ii,2}(index,1:2); % record value that produces optimal MI for temperature
% total{ii,10}(1,4) = index;
% [total{ii,10}(2,1) , index] = max (total{ii,5}(:,7));  % MI for 57 C isotherm
% total{ii,10}(2,2:3) = total{ii,2}(index,1:2); % record value that produces optimal MI for 57 C isotherm
% total{ii,10}(2,4) = index;
% 
% % Record optimal number of false pixels for 57 C isotherm
% total{ii,11} = zeros(1,4);
% [total{ii,11}(1,1) , index] = min (total{ii,6}(:,7,3));  % False pixel number for 57 C isotherm
% total{ii,11}(1,2:3) = total{ii,2}(index,1:2); % record value that produces number of false pixels
% total{ii,11}(1,4) = index;

end
        
%         % Record optimal L information
%         [total{ii,7}(1,1) , index] = min (total{ii,2}(:,2));  % record optimal L1
%         total{ii,7}(1,2) = total{ii,2}(index,1);   % record value that produces optimal L1
%         total{ii,7}(1,3) = index; % record index that produces optimal L1
%         [total{ii,7}(2,1) , index] = min (total{ii,2}(:,3));  % record optimal L2
%         total{ii,7}(2,2) = total{ii,2}(index,1);   % record value that produces optimal L2
%         total{ii,7}(2,3) = index; % record index that produces optimal L2
%         [total{ii,7}(3,1) , index] = min (total{ii,2}(:,4));  % record optimal L_inf
%         total{ii,7}(3,2) = total{ii,2}(index,1);   % record value that produces optimal L_inf
%         total{ii,7}(3,3) = index; % record index that produces optimal L_inf
%         
%         % Record optimal 57 C isotherm DSC information
%         [total{ii,8}(1) , index] = max (total{ii,3}(:,7));  % record optimal Dice
%         total{ii,8}(2) = total{ii,2}(index,1);   % record value that produces optimal Dice
%         total{ii,8}(3) = index; % record index that produces optimal Dice
%         
%         % Record optimal 57 C Hausdorff distance information
%         [total{ii,9}(1) , index] = min (total{ii,4}(:,7)); % record optimal Hausdorff Distance
%         total{ii,9}(2) = total{ii,2}(index,1); % record value that produces optimal Hausdorff distance
%         total{ii,9}(3) = index;
%         
%         % Record optimal temperature MI and 57 C isotherm MI
%         [total{ii,10}(1,1) , index] = max (total{ii,2}(:,5));  % MI for temperature
%         total{ii,10}(1,2) = total{ii,2}(index,1); % record value that produces optimal MI for temperature
%         total{ii,10}(1,3) = index;
%         [total{ii,10}(2,1) , index] = max (total{ii,5}(:,7));  % MI for 57 C isotherm
%         total{ii,10}(2,2) = total{ii,2}(index,1); % record value that produces optimal MI for 57 C isotherm
%         total{ii,10}(2,3) = index;
%         
%         % Record optimal number of false pixels for 57 C isotherm
%         [total{ii,11}(1,1) , index] = min (total{ii,6}(:,7,3));  % False pixel number for 57 C isotherm
%         total{ii,11}(1,2) = total{ii,2}(index,1); % record value that produces number of false pixels
%         total{ii,11}(1,3) = index;
%     end
end


num_runs = size(total_all,1);
num_var   = size(total_all{1,2},1);
L1 = zeros( num_var, num_runs);
L2 = L1;
L_inf = L1;
dice = L1;
HD = L1;
MI_temp = L1;
MI_57 = L1;
false_pix = L1;

for ii = 1:num_runs
    
    L1(:,ii)        = total_all{ii,2}(:,2);
    L2(:,ii)        = total_all{ii,2}(:,3);
    L_inf(:,ii)     = total_all{ii,2}(:,4);
    dice(:,ii)      = total_all{ii,3}(:,7);
    HD(:,ii)        = total_all{ii,4}(:,7);
    MI_temp(:,ii)   = total_all{ii,2}(:,5);
    MI_57(:,ii)     = total_all{ii,5}(:,7);
    false_pix(:,ii) = total_all{ii,6}(:,7,3);
    
end
clear ii

end
% % Descritpive stats of all...
% L1_mean = mean( L1,2);   % Plot all the L norms together; plot
% L1_median = median( L1,2);
% L1_stDev = std(L1')';
% L2_mean = mean( L2,2);
% L2_median = median( L2,2);
% L2_stDev  = std ( L2')';
% L_inf_mean = mean( L_inf,2);
% L_inf_median = median( L_inf,2);
% L_inf_stDev  = std ( L_inf')';
% dice_mean = mean( dice,2);     % Plot DSC and HD together; plotyy
% dice_median = median( dice,2);
% dice_stDev  = std ( dice')';
% HD_mean = mean( HD,2);
% HD_median = median( HD,2);
% HD_stDev  = std ( HD')';
% MI_temp_mean = mean( MI_temp,2);       % Plot MI together;
% MI_temp_median = median( MI_temp,2);
% MI_temp_stDev  = std ( MI_temp')';
% MI_57_mean = mean( MI_57,2);
% MI_57_median = median( MI_57,2);
% MI_57_stDev  = std ( MI_57')';
% false_pix_mean = mean( false_pix,2);
% false_pix_median = median( false_pix,2);
% false_pix_stDev  = std ( false_pix')';
