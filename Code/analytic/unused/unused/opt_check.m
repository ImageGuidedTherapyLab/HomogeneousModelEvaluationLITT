% This script finds the best mu_eff for the different studies.
clear
tic;
% cell_data = csvimport('datasummary.txt');
% headers = cell_data(1,1:3);
% mu_eff_data = cell2mat(cell_data(2:end,:));

% Identify the studies to be examined.
opttype = 'bestfit49' ;

Study_paths = cell (1,2);
% Study_paths {1,1} = 'Study0035';
% Study_paths {1,2} = '0530';
% Study_paths {2,1} = 'Study0030';
% Study_paths {2,2} = '0495';
% Study_paths {3,1} = 'Study0030';
% Study_paths {3,2} = '0497';
% Study_paths {4,1} = 'Study0030';
% Study_paths {4,2} = '0491';
% Study_paths {5,1} = 'Study0030';
% Study_paths {5,2} = '0496';
% Study_paths {6,1} = 'Study0030';
% Study_paths {6,2} = '0490';
% Study_paths {7,1} = 'Study0017';
% Study_paths {7,2} = '0378';
% Study_paths {8,1} = 'Study0025';
% Study_paths {8,2} = '0438';
% Study_paths {9,1} = 'Study0025';
% Study_paths {9,2} = '0435';
% Study_paths {10,1} = 'Study0025';
% Study_paths {10,2} = '0440';
% Study_paths {11,1} = 'Study0025';
% Study_paths {11,2} = '0436';
% Study_paths {12,1} = 'Study0028';
% Study_paths {12,2} = '0466';
% Study_paths {13,1} = 'Study0028';
% Study_paths {13,2} = '0468';
% Study_paths {14,1} = 'Study0028';
% Study_paths {14,2} = '0471';
% Study_paths {15,1} = 'Study0026';
% Study_paths {15,2} = '0447';
% Study_paths {16,1} = 'Study0026';
% Study_paths {16,2} = '0457';
Study_paths {1,1} = 'Study0030';
Study_paths {1,2} = '0497';
Study_paths {2,1} = 'Study0018';
Study_paths {2,2} = '0402';
Study_paths {3,1} = 'Study0029';
Study_paths {3,2} = '0476';
Study_paths {4,1} = 'Study0025';
Study_paths {4,2} = '0436';
Study_paths {5,1} = 'Study0028';
Study_paths {5,2} = '0466';
Study_paths {6,1} = 'Study0028';
Study_paths {6,2} = '0468';
Study_paths {7,1} = 'Study0028';
Study_paths {7,2} = '0471';
Study_paths {8,1} = 'Study0021';
Study_paths {8,2} = '0409';

datasummary = dlmread('datasummaryL2_10sourceNewton49.txt',',',1,0);

num_total = size(datasummary,1);
num_studies = size(Study_paths,1);
best_iter = zeros(num_studies,1);
index = zeros(num_studies,1);
for jj = 1:num_studies
    
    for ii = 1:num_total
        
        TF = strcmp(strcat('0',num2str(datasummary(ii,2))),Study_paths{jj,2});
        
        if TF ==1
            index(jj)=ii;
            
        end
    end
    best_iter(jj) = datasummary(index(jj),3);

end

clear ii jj

[L2norm, dice, tmap_model, MRTI_crop,mat_struct,quality] = check_opt_GoodBadUgly( Study_paths, 'bestfit49', best_iter);

% 
% for ii = 1:num_studies
%     disp('Start ')
%     disp(num2str(ii))
%     [L2norm, dice, tmap_model, MRTI_crop,mat_struct,quality] = check_opt_GoodBadUgly( Study_paths, 'bestfit49', datasummary(:,3));
%     
%     total_path{ii,1} = strcat(Study_paths{ii,1}, '/', Study_paths{ii,2});
%     input_path{1,1} = Study_paths{ii,1};
%     input_path{1,2} = Study_paths{ii,2};
%     %[ total_path{ii,2}, total_path{ii,3} ] = Check_ablation44 ( input_path , opttype);
%     %[ total_path{ii,2}, total_path{ii,3} ] = Check_ablation55 ( input_path , opttype);
%     %[ total_path{ii,2}, total_path{ii,3} ] = Check_ablation66 ( input_path , opttype);
%     [ total_path{ii,2}, total_path{ii,3} ] = Check_ablation77 ( input_path , opttype);
%     %[ total_path{ii,2}, total_path{ii,3} , model_temp_dice_shoulder, MRTI_crop_dice_shoulder] = Check_ablation55 ( input_path , opttype);
% end
%[ H0, H1, dice_values ] = Check_ablation ( Study_paths, mu_eff_opt );
toc

% model_deg_threshold = model_temp_dice_shoulder >= 57;
% MRTI_deg_threshold = MRTI_crop_dice_shoulder >= 57;
% n_model = sum(sum( model_deg_threshold ));
% n_MRTI = sum(sum( MRTI_deg_threshold ));
% 
% MRTI_label = MRTI_deg_threshold .* 2;
% 
% intersection = model_deg_threshold + MRTI_deg_threshold;
% intersection = intersection > 1;
% intersection_label = model_deg_threshold + MRTI_label;
% n_intersection = sum(sum( intersection ));
% dice = 2*n_intersection / (n_model + n_MRTI);
% figure(4); imagesc(model_deg_threshold);
% figure(5); imagesc(MRTI_deg_threshold);
% figure(6); imagesc(intersection_label, [ 0 3]);
% figure(7); imagesc(model_temp_dice_shoulder,[30 90]);
% figure(8); imagesc(MRTI_crop_dice_shoulder, [30 90]);