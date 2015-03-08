% This script finds the best mu_eff for the different studies.
function metric_generation (choice,quick_choice);
%choice = 1; % 1 = mu; 2 = perf; 3 = cond;
tic

cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file

opttype = 'bestfit50';
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);
for ii = 1:num_studies
    
    Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end
clear ii


cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/libraries

if choice == 1
    
    %load ('all_opt_mu.mat' );
    load ('all_opt_mu_small.mat' );
    
elseif choice == 2
    
    load ('all_opt_perf.mat' );
    
elseif choice == 3
    
    load ('all_opt_cond.mat' );
    
elseif choice == 4
    
    %load ('all_opt_perf_mu_400_short.mat' );
    load ('opt_perf_mu_400_long.mat');
    
elseif choice == 5
    
    load ('rand_opt_perf_mu.mat' );
    
end
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
% indexC = strfind(opt.paths,'Study0035');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom,:)=[];
% datasummary(toss_index_phantom,:)=[];
% num_studies = size(datasummary,1);
%
% indexC = strfind(opt.paths,'0457');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom-num_studies,:)=[];
% datasummary(toss_index_phantom-num_studies,:)=[];
% num_studies = size(datasummary,1);
%
% indexC = strfind(opt.paths,'0476');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom-num_studies,:)=[];
% datasummary(toss_index_phantom-num_studies,:)=[];
% num_studies = size(datasummary,1);
%
% indexC = strfind(opt.paths,'0436');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom-num_studies,:)=[];
% datasummary(toss_index_phantom-num_studies,:)=[];

% 0491 is suspect: small ROI, ablation fills entire ROI (ie not good)
% 0490 is suspect: small ROI, ablation fills entire ROI
% 0378 is suspect: has very funny shape
% 0385 is really small (but nice and round)
% 0476 is really small (and funny shaped)
% 0435 is really funny shaped
% 0436 small
% 0466 small
% 0468 Really small
% 0471 small and funny shaped
% 0477 small
% 0450 almost small

% From mu_eff_data, find the matching study's(ies') mu_eff value(s)
num_studies = size(Study_paths,1);

% clear ii
total = cell(num_studies+1,11);

% Labels for the sundry cell columns
total{1,1}    = 'Study';
total{1,2}{1} = '1= Parameter value';
total{1,2}{2} = '2:4 = L1, L2, L_inf norms';
total{1,2}{3} = '5:7 = Temp MI, Total heating, Max temp';
total{1,3}    = 'DSC isotherms 51:65; column index 7 is 57';
total{1,4}    = 'HD for isotherms';
total{1,5}    = 'MI for isotherms';
total{1,6}    = 'False pixel count for isotherms';
total{1,7}    = 'Optimal L norms;';
total{1,8}    = 'Optimal DSC for 57 C';
total{1,9}    = 'Optimal HD for 57 C';
total{1,10}   = 'Optimal temp and 57 C MI';
total{1,11}   = 'Optimal false pixel count for 57 C';

for ii = 2:(num_studies+1)
    % Display run information
    disp('Start ')
    disp(strcat (num2str(ii-1),' of ', num2str(num_studies)))
    total{ii,1} = strcat(Study_paths{ii-1,1}, '/', Study_paths{ii-1,2});
    fprintf('iter %s \n', total{ii,1});
    toc
    
    % Load and organize MRTI
    path_base = strcat ( 'workdir/',Study_paths{ii-1,1}, '/', Study_paths{ii-1,2}, '/opt');
    load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );
    [MRTI_crop, MRTI_isotherm, MRTI_list,n_MRTI] = organize_MRTI( inputdatavars );
    
    % Generate model data
    [model_crop] = organize_model( inputdatavars, all_opt_fig, no_pwr_fig,sim_dim, choice );
    
    if choice == 5|| choice==4
        [ total{ii,2}, total{ii,3}, total{ii,4}, total{ii,5}, total{ii,6} ] = metric_calculator_rand ( MRTI_crop, MRTI_isotherm, MRTI_list,n_MRTI,model_crop,inputdatavars,summary,quick_choice );
        [total] = optimal_metrics_rand(total,ii,quick_choice);
    else        
        [ total{ii,2}, total{ii,3}, total{ii,4}, total{ii,5}, total{ii,6} ] = metric_calculator ( MRTI_crop, MRTI_isotherm, MRTI_list,n_MRTI,model_crop,inputdatavars,summary );
        [total] = optimal_metrics (total,ii);
    end
    

%     % Record optimal L information
%     total{ii,7} = zeros(3);
%     [total{ii,7}(1,1) , index] = min (total{ii,2}(:,2));  % record optimal L1
%     total{ii,7}(1,2) = total{ii,2}(index,1);   % record value that produces optimal L1
%     total{ii,7}(1,3) = index; % record index that produces optimal L1
%     [total{ii,7}(2,1) , index] = min (total{ii,2}(:,3));  % record optimal L2
%     total{ii,7}(2,2) = total{ii,2}(index,1);   % record value that produces optimal L2
%     total{ii,7}(2,3) = index; % record index that produces optimal L2
%     [total{ii,7}(3,1) , index] = min (total{ii,2}(:,4));  % record optimal L_inf
%     total{ii,7}(3,2) = total{ii,2}(index,1);   % record value that produces optimal L_inf
%     total{ii,7}(3,3) = index; % record index that produces optimal L_inf
%     
%     % Record optimal 57 C isotherm DSC information
%     total{ii,8} = zeros(1,3);
%     [total{ii,8}(1) , index] = max (total{ii,3}(:,7));  % record optimal Dice
%     total{ii,8}(2) = total{ii,2}(index,1);   % record value that produces optimal Dice
%     total{ii,8}(3) = index; % record index that produces optimal Dice
%     
%     % Record optimal 57 C Hausdorff distance information
%     total{ii,9} = zeros(1,3);
%     [total{ii,9}(1) , index] = min (total{ii,4}(:,7)); % record optimal Hausdorff Distance
%     total{ii,9}(2) = total{ii,2}(index,1); % record value that produces optimal Hausdorff distance
%     total{ii,9}(3) = index;
%     
%     % Record optimal temperature MI and 57 C isotherm MI
%     total{ii,10} = zeros(2,3);
%     [total{ii,10}(1,1) , index] = max (total{ii,2}(:,5));  % MI for temperature
%     total{ii,10}(1,2) = total{ii,2}(index,1); % record value that produces optimal MI for temperature
%     total{ii,10}(1,3) = index;
%     [total{ii,10}(2,1) , index] = max (total{ii,5}(:,7));  % MI for 57 C isotherm
%     total{ii,10}(2,2) = total{ii,2}(index,1); % record value that produces optimal MI for 57 C isotherm
%     total{ii,10}(2,3) = index;
%     
%     % Record optimal number of false pixels for 57 C isotherm
%     total{ii,11} = zeros(1,3);
%     [total{ii,11}(1,1) , index] = min (total{ii,6}(:,7,3));  % False pixel number for 57 C isotherm
%     total{ii,11}(1,2) = total{ii,2}(index,1); % record value that produces number of false pixels
%     total{ii,11}(1,3) = index;
    
end
%[ H0, H1, dice_values ] = Check_ablation ( Study_paths, mu_eff_opt );
toc

% Save
cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/libraries

if choice == 1
    
    save ('GPU_dict_mu.mat','total','summary','-v7.3');
    
elseif choice == 2
    
    save ('GPU_dict_perf.mat','total','summary','-v7.3');
    
elseif choice == 3
    
    save ('GPU_dict_cond.mat','total','summary','-v7.3');
    
elseif choice == 4
    
    save ('GPU_dict_perf_mu_global_400','total','summary','-v7.3');
    
elseif choice == 5
    
    save ('GPU_dict_perf_mu_rand','total','summary','-v7.3');
    
end
 
end