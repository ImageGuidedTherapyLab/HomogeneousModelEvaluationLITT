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

% This silly script is a sensitivity study. I put some of the data in
% /FUS4/data2/sjfahrenholtz/Matlab/Tests

function [ total, dice,hd, summary ] = Check_ablation_choice ( Study_paths, opttype, choice );

% Make the LOOCV iteration system
n_patients = size( Study_paths,1); % This is the number of patients

threshold_temps = 51:65;
num_threshold_temps = length(threshold_temps);


if choice ==1
    
    mu_eff(1) = 0.008;
    mu_eff(2:101) = linspace(101,200,100);
    %mu_eff(2:5001) = linspace(1,5000,5000);
    
    dice = zeros( length(mu_eff),num_threshold_temps); % Initialize the number of DSC (dice) values
    
elseif choice ==2
    
    w_perf = 0.01;
    %w_perf(2:35) = linspace ( 0.5, 17, 34);
    %w_perf(2:34) = linspace ( 0.5, 16.5,33);
    w_perf(2:133) = linspace ( 0.125, 16.5,132);
    
    dice = zeros( length(w_perf),num_threshold_temps); % Initialize the number of DSC (dice) values
    
elseif choice ==3
    
    %k_cond = linspace ( 0.01,2,200);
    k_cond = linspace ( 0.52, 3, 249);
    
    dice = zeros( length(k_cond),num_threshold_temps); % Initialize the number of DSC (dice) values
    
end

path_base = strcat ( 'workdir/',Study_paths{1,1}, '/', Study_paths{1,2}, '/opt');
load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );
mu = str2num( inputdatavars.cv.mu_eff_healthy );
kk = inputdatavars.cv.k_0;
ww = inputdatavars.cv.w_0;
%for ii = 1
% This section writes a new TmpDataInput.mat file for later reading.
for ii = 1:n_patients
    
    if choice ==1   % Mu
        %k_cond = 0.527;
        %w_perf = 6;
        k_cond = kk;
        w_perf = ww;
        summary.k_cond = k_cond;
        summary.w_perf = w_perf;
        [total,dice, hd] = temperature_obj_fxn_GPU_choice ( inputdatavars, 25, mu_eff, w_perf, k_cond, choice );
        
    elseif choice ==2   % perf
        % mu_eff = 180;
        % k_cond = 0.527;
        mu_eff = 180;
        k_cond = kk;
        summary.mu_eff = mu_eff;
        summary.k_cond = k_cond;
        [total,dice, hd] = temperature_obj_fxn_GPU_choice ( inputdatavars, 25, mu_eff, w_perf, k_cond, choice );
        
    elseif choice ==3    % cond
        % mu_eff = 180;
        % w_perf = 6;
        mu_eff = 180;
        w_perf = ww;
        summary.mu_eff = mu_eff;
        summary.w_perf = w_perf;
        [total,dice, hd] = temperature_obj_fxn_GPU_choice ( inputdatavars, 25, mu_eff, w_perf, k_cond, choice );
    end
  
end

end