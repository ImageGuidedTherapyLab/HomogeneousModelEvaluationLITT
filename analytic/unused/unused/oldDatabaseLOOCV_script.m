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

Patient_paths     = cell (5,2);
Patient_paths {1,1} = 'Patient0002';
Patient_paths {1,2} = '000';
Patient_paths {2,1} = 'Patient0002';
Patient_paths {2,2} = '001';
Patient_paths {3,1} = 'Patient0002';
Patient_paths {3,2} = '021';
Patient_paths {4,1} = 'Patient0006';
Patient_paths {4,2} = '007';
Patient_paths {5,1} = 'Patient0006';
Patient_paths {5,2} = '009';

mu_eff_opt = zeros (5,1); % Initialize a couple looped variables
min_obj_fxn_indices = zeros (5,1);
% for ii = 1:5 % This loop finds the optimized mu_eff values from datasets that have already been optimized
for ii = 1:5
    single_path = strcat( 'workdir/', Patient_paths{ii,1}, '/', Patient_paths{ii,2}, '/opt/');
    tab_dat = load (strcat( single_path,'dakota_q_newton_heating.in.tabular.dat') );
    [~,min_obj_fxn_indices(ii)] = min( tab_dat(:,3) );
    mu_eff_opt (ii) = tab_dat( min_obj_fxn_indices(ii),2);
end

[ H0, H1 ] = LOOCV_t_test ( Patient_paths, mu_eff_opt );

H0
H1