% This script finds the best mu_eff for the different studies.
cell_data = csvimport('datasummary.txt');
headers = cell_data(1,1:3);
mu_eff_data = cell2mat(cell_data(2:end,:));

% Identify the studies to be examined.
Study_paths = cell (1,2);
Study_paths {1,1} = 'Study0035';
Study_paths {1,2} = '0530';

% From mu_eff_data, find the matching study's(ies') mu_eff value(s)
num_studies = size(Study_paths,1);
matching_num = zeros(1,num_studies);
mu_eff_index = zeros(1,num_studies);
mu_eff_opt   = zeros(1,num_studies);
for ii = 1:num_studies
    
    matching_num(ii) = str2num(Study_paths{(ii),2});
    mu_eff_index(ii) = find( mu_eff_data(:,1) == matching_num(ii));
    mu_eff_opt(ii) = mu_eff_data(2,mu_eff_index(ii));
end

[ H0, H1 ] = Check_ablation ( Study_paths, mu_eff_opt );