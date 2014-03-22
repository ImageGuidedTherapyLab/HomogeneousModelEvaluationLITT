% I have the list of good runs, now I need to code it. 
% 
% The input: vector of paths to whatever interested datasets and vector of optimized mu_eff values.
% (I assume that the inverse problem is run before the LOOCV begins).
% 
% The process:
% Its iterative. Heres one iteration of leave-one-out cross validation (LOOCV). Average the optimized mu_eff for all but one ablation. 
% Use the averaged optimized mu_eff for the reserved ablation. 
% That means a dakota.in file will be written with the default constants and the averaged optimized mu_eff value.
% That dakota.in file wile be read by the BHTE code, the BHTE code runs, and a prediction is made. 
% Then, the prediction 57 deg_C isotherm is compared to the MRTI 57 dec_C isotherm using a Dice similarity coefficient (DSC). That is one iteration.
% Overall process iterates the LOOCV algorithm such that there are n possible iterations for n ablation datasets resulting in n DSC values. 
% The DSC values mean and variance are calculated. The mean and variance are sent to the t-test to find the hypotheses results.
% 
% The output: The output is a binary acceptance/rejection of the null and alternative hypotheses.

% This script finds the best mu_eff for the different studies.
MATLAB_on_off = 1; % 0 = off/false; 1 = on/true
opt_type = 'heating' ;
cell_data = csvimport('datasummary.txt');
headers = cell_data(1,1:3);
mu_eff_data = cell2mat(cell_data(2:end,:));

% Identify the studies to be examined.
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

num_studies = size(Study_paths,1);
matching_num = zeros(1,num_studies);
mu_eff_index = zeros(1,num_studies);
mu_eff_opt   = zeros(1,num_studies);
for ii = 1:num_studies
    matching_num(ii) = str2num(Study_paths{(ii),2});
    mu_eff_index(ii) = find( mu_eff_data(:,1) == matching_num(ii));
    mu_eff_opt(ii) = mu_eff_data(mu_eff_index(ii),2);
end

if MATLAB_on_off == 0 % Make sure brainsearch.py has MATLAB flag off ('false')
    [ hh, dice_values] = LOOCV_t_test_DF ( Study_paths, mu_eff_opt, opt_type );
    
elseif MATLAB_on_off == 1 % Make sure brainsearch.py has MATLAB flag on ('true')
    [ hh, dice_values ] = LOOCV_t_test ( Study_paths, mu_eff_opt, opt_type ); 
    
else
    hh = 'ERROR: MATLAB on/off flag is not 0 or 1.' ;
    
end

hh