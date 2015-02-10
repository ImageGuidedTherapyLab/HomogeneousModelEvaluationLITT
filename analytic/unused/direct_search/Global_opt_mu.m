% This script finds the best mu_eff for the different studies.
tic;
clear
close all
clc

% Identify the studies to be examined.
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file

opttype = 'bestfit50' ;
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);

opt.paths = cell (num_studies,2);
for ii = 1:num_studies
    
    Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end

clear ii
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

% From mu_eff_data, find the matching study's(ies') mu_eff value(s)
num_studies = size(Study_paths,1);
% matching_num = zeros(1,num_studies);
% mu_eff_index = zeros(1,num_studies);
% mu_eff_opt   = zeros(1,num_studies);
% for ii = 1:num_studies
%     
%     matching_num(ii) = str2num(Study_paths{(ii),2});
%     mu_eff_index(ii) = find( mu_eff_data(:,1) == matching_num(ii));
%     mu_eff_opt(ii) = mu_eff_data(mu_eff_index(ii),2);
% end
% clear ii
total = cell(num_studies,7);
input_path = cell(1,2);

for ii = 1:num_studies
    % Display run information
    disp('Start ')
    disp(strcat (num2str(ii),' of ', num2str(num_studies)))
    total{ii,1} = strcat(Study_paths{ii,1}, '/', Study_paths{ii,2});
    fprintf('iter %s \n', total{ii,1});
    
    % Do global optimization
    input_path{1,1} = Study_paths{ii,1};
    input_path{1,2} = Study_paths{ii,2};
    [ total{ii,2}, total{ii,3}, total{ii,4} ] = Check_ablation_mu ( input_path , opttype);
    
    % Record optimal L2 information
    total{ii,5} = zeros(1,3);
    [total{ii,5}(1) , index] = min (total{ii,2}(:,2));  % record optimal L2
    total{ii,5}(2) = total{ii,2}(index,1);   % record mu_eff that produces optimal L2
    total{ii,5}(3) = index; % record index that produces optimal L2
    
    % Record optimal 57 C isotherm DSC information
    total{ii,6} = zeros(1,3);
    [total{ii,6}(1) , index] = max (total{ii,3}(:,7));  % record optimal Dice
    total{ii,6}(2) = total{ii,2}(index,1);   % record mu_eff that produces optimal Dice
    total{ii,6}(3) = index; % record index that produces optimal Dice
    
    % Record optimal 57 C Hausdorff distance information
    total{ii,7} = zeros(1,3);
    [total{ii,7}(1) , index] = min (total{ii,4}(:,7)); % record optimal Hausdorff Distance
    total{ii,7}(2) = total{ii,2}(index,1); % record mu_eff that produces optimal Hausdorff distance
    total{ii,7}(3) = index;
end
%[ H0, H1, dice_values ] = Check_ablation ( Study_paths, mu_eff_opt );
toc

% Save
cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search
save ('GPU_global_mu.mat','total')

cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

data_filename = 'datasummaryL2_10sourceNewton50.txt';
data_filename_out= 'data_summary_GPU_hd.txt';

fin  = fopen(data_filename);
fout = fopen( data_filename_out,'w');



%% This records the optimal Dice
f_line = fgetl(fin); % basically skip first line
fprintf(fout,'%s\n',f_line);
fclose(fin);
for ii = 1:num_studies
    
    first_half  = strcat( Study_paths{ii,1}(6:end),',',Study_paths{ii,2},',1,', num2str( total{ii,6}(2)) , ',1.3854600e-07,','0.00e+00,'); % beginning to robin
    second_half = strcat( num2str( total{ii,6}(1)),',', num2str( total{ii,5}(1)),',', num2str( total{ii,5}(1)),',', num2str( 1-total{ii,6}(1) )); % obj fxns
    f_line = strcat( first_half,second_half);    
    fprintf( fout, '%s\n', f_line);

end

fclose(fout);



data_filename_out= 'data_summary_GPU_L2.txt';

fin  = fopen(data_filename);
fout = fopen( data_filename_out,'w');



%% This records the optimal L2
f_line = fgetl(fin); % basically skip first line
fprintf(fout,'%s\n',f_line);
fclose(fin);
for ii = 1:num_studies
    
    first_half  = strcat( Study_paths{ii,1}(6:end),',',Study_paths{ii,2},',1,', num2str( total{ii,4}(2)) , ',1.3854600e-07,','0.00e+00,'); % beginning to robin
    second_half = strcat( num2str( total{ii,3}( total{ii,4}(2), 7)),',', num2str( total{ii,4}(1)),',', num2str( total{ii,4}(1)),',', num2str( 1-total{ii,3}( total{ii,4}(2), 7) )); % obj fxns
    f_line = strcat( first_half,second_half);    
    fprintf( fout, '%s\n', f_line);

end

fclose(fout);
