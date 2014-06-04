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

% H0 is null hypothesis
% H1 is alternative hypothesis
% Study_paths is the paths to the ablations that are being run through
%   the LOOCV algorithm
% mu_eff_opt is a vector of the inverse problem optimized mu_eff values in
%   1/m units

function [ hh, dice_values ] = LOOCV_t_test_DF ( Study_paths, mu_eff_opt, alpha_opt, opt_type );

% Make the LOOCV iteration system
n_patients = length( mu_eff_opt); % This is the number of patients
% n_patients = 1;
dice_values = zeros( n_patients,1); % Initialize the number of DSC (dice) values
for ii = 1:n_patients
    % Set up LOOCV for mu_eff
    mu_eff_iter = mu_eff_opt; % Make a copy of both the mu_eff values and the paths
    mu_eff_iter ( ii ) = []; % Remove
    mu_eff_iter = mean ( mu_eff_iter );           % Add alpha to LOOCV here; write *.in.* file; run brainsearch.py; brainsearch.py makes *.out.* file; read in *.out.* file
%     mu_s_p = params_iter.cv.mu_s * ( 1 - params_iter.cv.anfact );
%     mu_a_iter = (-3*mu_s_p + sqrt( 9*mu_s_p^2 + 12 * mu_eff_iter^2))/6; % Used quadratic equation to solve for mu_a in terms of g, mu_s, and mu_eff
    
    alpha_iter = alpha_opt;
    alpha_iter (ii) = [];
    alpha_iter = mean (alpha_iter);
    
    input_file = cell (18,1);
    input_file {1,1} = '12 variables';
    input_file {2,1} = strcat ( num2str(mu_eff_iter),' mu_eff_healthy');
    input_file {3,1} = strcat ( num2str(alpha_iter), ' alpha_healthy');
    input_file {4,1} = '-9.000000000000000e+01 y_rotate';
    input_file {5,1} = '2.287330000000000e-05 gamma_healthy';
    input_file {6,1} = '1.600000000000000e-01 x_displace';
    input_file {7,1} = '1.450000000000000e-01 y_displace';
    input_file {8,1} = '0.000000000000000e+00 z_displace';
    input_file {9,1} = '3.700000000000000e+01 body_temp';
    input_file{10,1} = '0.000000000000000e+00 x_rotate';
    input_file{11,1} = '0.000000000000000e+00 z_rotate';
    input_file{12,1} = '2.100000000000000e+01 probe_init';
    input_file{13,1} = '0.000000000000000e+00 robin_coeff';
    input_file{14,1} = '1 functions';
    input_file{15,1} = '1 ASV_1:obj_fn';
    input_file{16,1} = '1 derivative_variables';
    input_file{17,1} = '1 DVV_1:mu_eff_healthy';
    input_file{18,1} = '0 analysis_components';
    
    
    
    
    
    
    % This section prepares the varied parameters into a .mat file for the
    % thermal code to run
    param_file = strcat( 'workdir/', Study_paths { ii,1 }, '/', Study_paths { ii,2 }, '/opt/', opt_type, '.in.1');
    python_command = strcat( 'unix(''python ./brainsearch.py --param_file ./', param_file, ''')');   % unix(''python test_saveFile.py'')
    evalc(python_command); 
   
    % This section runs the thermal code
%     [metric, dice_iter, thermal_model, MRTI_crop] = fast_temperature_obj_fxn_sanity ( params_iter, 1 );
    dice_values (ii) = dice_iter ;
    clear mu_eff_iter;
end

[hh.H, hh.ptest] = ttest( dice_values, 0.7, 0.05, 'right');

% mean_dice = mean ( dice_values ); % Average Dice values
% std_dice = std ( dice_values ); % Stardard deviation of Dice values.
%test_statistic = ( mean_dice - 0.7 ) ./ ( std_dice ./ sqrt ( n_patients )); % Calculate the test_statistic
%p_value = tcdf ( test_statistic , n_patients , 'lower' );
% H0.p_value = p_value;
% H0.result = 1;
% if p_value <= 0.05
%     H0.result = 0 ;
% else
%     H0.result = 1;
% end

end