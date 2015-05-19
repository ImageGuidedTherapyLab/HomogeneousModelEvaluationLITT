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

function [ hh, dice_values ] = naive_t_test_DiceTemp ( Study_paths, mu_eff_opt, alpha_opt, best_iter, opttype);

% Make the LOOCV iteration system
n_patients = length( mu_eff_opt); % This is the number of patients
% n_patients = 1;
%path_iter = cell(1,2);
L2norm = zeros( n_patients,1);
dice_values = zeros( n_patients,1); % Initialize the number of DSC (dice) values
for ii = 1:n_patients
    
    disp( fprintf('%d of %d\n',ii,n_patients) );
    
    % Set up LOOCV for mu_eff
    mu_eff_iter = mu_eff_opt; % Make a copy of both the mu_eff values and the paths
    mu_eff_iter ( ii ) = []; % Remove the 0
    mu_eff_iter = 180;           % Add alpha to LOOCV here; write *.in.* file; run brainsearch.py; brainsearch.py makes *.out.* file; read in *.out.* file
    %mu_s_p = params_iter.cv.mu_s * ( 1 - params_iter.cv.anfact );
    %mu_a_iter = (-3*mu_s_p + sqrt( 9*mu_s_p^2 + 12 * mu_eff_iter^2))/6; % Used quadratic equation to solve for mu_a in terms of g, mu_s, and mu_eff
    
    % Set up LOOCV for alpha
    alpha_iter = alpha_opt;
    alpha_iter (ii) = [];
    alpha_iter = mean ( alpha_iter );
    
    DAKOTA_in_writer_naive ( Study_paths(ii,:), mu_eff_iter, alpha_iter, best_iter(ii), opttype );
    % This section prepares the varied parameters into a .mat file for the
    % thermal code to run
    %     param_file  = strcat( 'workdir/', Study_paths { ii,1 }, '/', Study_paths { ii,2 }, '/opt/optpp_pds.', opt_type, '.in.1');
    %     result_file = strcat( 'workdir/', Study_paths { ii,1 }, '/', Study_paths { ii,2 }, '/opt/optpp_pds.', opt_type, '.out.1');
    path_base = strcat( 'workdir/',Study_paths{ii,1},'/',Study_paths{ii,2},'/opt/');
    param_file = strcat( path_base,'optpp_pds.naive.in.1');
    result_file = strcat( path_base, 'optpp_pds.naive.out.1');
    python_command = strcat( '[errorcode, codestdout] = unix(''python ./brainsearch.py --param_file ./', param_file, ' ./', result_file, ''')');   % unix(''python test_saveFile.py'')
    disp(python_command);
    evalc(python_command)
    
    output = dlmread( strcat( './',result_file));
    L2norm (ii) = output(1);
    dice_values (ii) = output(3);

end

[hh.H, hh.ptest, hh.ci, hh.stats] = ttest( dice_values, 0.7, 0.05, 'right');

% mean_dice = mean ( dice_values ); % Average Dice values
% std_dice = std ( dice_values ); % Stardard deviation of Dice values.
% test_statistic = ( mean_dice - 0.7 ) ./ ( std_dice ./ sqrt ( n_patients )); % Calculate the test_statistic
% p_value = tcdf ( test_statistic , n_patients , 'lower' );
% H0.p_value = p_value;
% H0.result = 1;
% if p_value <= 0.05
%     H0.result = 0 ;
% else
%     H0.result = 1;
% end

end
