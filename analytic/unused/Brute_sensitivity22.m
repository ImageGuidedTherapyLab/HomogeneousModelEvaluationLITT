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
Study_paths {1,1} = 'Study0023';
Study_paths {1,2} = '0428';
% Study_paths {18,1} = 'Study0026';
% Study_paths {18,2} = '0453';
% Study_paths {19,1} = 'Study0026';
% Study_paths {19,2} = '0450';
% Study_paths {20,1} = 'Study0026';
% Study_paths {20,2} = '0451';
% Study_paths {21,1} = 'Study0022';
% Study_paths {21,2} = '0418';
% Study_paths {22,1} = 'Study0022';
% Study_paths {22,2} = '0417';
% Study_paths {23,1} = 'Study0021';
% Study_paths {23,2} = '0409';
% Study_paths {24,1} = 'Study0021';
% Study_paths {24,2} = '0414';
% Study_paths {25,1} = 'Study0021';
% Study_paths {25,2} = '0415';

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
total_path = cell(num_studies,3);
input_path = cell(1,2);

for ii = 1:num_studies
    disp('Start ')
    disp(num2str(ii))
    total_path{ii,1} = strcat(Study_paths{ii,1}, '/', Study_paths{ii,2});
    input_path{1,1} = Study_paths{ii,1};
    input_path{1,2} = Study_paths{ii,2};
    %[ total_path{ii,2}, total_path{ii,3} ] = Check_ablation44 ( input_path , opttype);
    %[ total_path{ii,2}, total_path{ii,3} ] = Check_ablation55 ( input_path , opttype);
    %[ total_path{ii,2}, total_path{ii,3} ] = Check_ablation66 ( input_path , opttype);
    [ total_path{ii,2}, total_path{ii,3} ] = Check_ablation77 ( input_path , opttype);
    %[ total_path{ii,2}, total_path{ii,3} , model_temp_dice_shoulder, MRTI_crop_dice_shoulder] = Check_ablation55 ( input_path , opttype);
end
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