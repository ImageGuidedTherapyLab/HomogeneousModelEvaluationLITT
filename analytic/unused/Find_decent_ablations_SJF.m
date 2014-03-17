% This script finds the best mu_eff for the different studies.
cell_data = csvimport('datasummary.txt');
headers = cell_data(1,1:3);
mu_eff_data = cell2mat(cell_data(2:end,:));

% Identify the studies to be examined.
Study_paths = cell (1,2);
Study_paths {1,1} = 'Study0035';
Study_paths {1,2} = '0530';
Study_paths {2,1} = 'Study0030';
Study_paths {2,2} = '0495';
Study_paths {3,1} = 'Study0030';
Study_paths {3,2} = '0497';
Study_paths {4,1} = 'Study0030';
Study_paths {4,2} = '0491';
Study_paths {5,1} = 'Study0030';
Study_paths {5,2} = '0496';
Study_paths {6,1} = 'Study0030';
Study_paths {6,2} = '0490';
Study_paths {7,1} = 'Study0017';
Study_paths {7,2} = '0378';
Study_paths {8,1} = 'Study0025';
Study_paths {8,2} = '0438';
Study_paths {9,1} = 'Study0025';
Study_paths {9,2} = '0435';
Study_paths {10,1} = 'Study0025';
Study_paths {10,2} = '0440';
% From mu_eff_data, find the matching study's(ies') mu_eff value(s)
num_studies = size(Study_paths,1);
matching_num = zeros(1,num_studies);
mu_eff_index = zeros(1,num_studies);
mu_eff_opt   = zeros(1,num_studies);
for ii = 1:num_studies
    
    matching_num(ii) = str2num(Study_paths{(ii),2});
    mu_eff_index(ii) = find( mu_eff_data(:,1) == matching_num(ii));
    mu_eff_opt(ii) = mu_eff_data(mu_eff_index(ii),2);
end

[ H0, H1 ] = Check_ablation ( Study_paths, mu_eff_opt );