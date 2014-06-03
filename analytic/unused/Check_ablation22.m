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

function [ total, dice_values ] = Check_ablation22 ( Study_paths, mu_eff_opt );

% Make the LOOCV iteration system
n_patients = size( Study_paths,1); % This is the number of patients
% n_patients = 1;
dice_values = zeros( n_patients,1); % Initialize the number of DSC (dice) values
for ii = 1
    % This section writes a new TmpDataInput.mat file for later reading.
    param_file = strcat( 'workdir/', Study_paths{ii,1}, '/', Study_paths{ii,2}, '/opt/optpp_pds.heating.in.1');
    python_command = strcat( 'unix(''python ./brainsearch.py --param_file ./', param_file, ''')');   % unix(''python test_saveFile.py'')
    evalc(python_command);
    
    % This section prepares the varied parameters into a .mat file for the
    % thermal code to run
    % Sample:  'python ./brainsearch.py --param_file '
    params_iter = load( 'TmpDataInput.mat' ); % Read in one dakota.in file to find the constant parameters
    params_iter.cv.mu_eff_healthy = mu_eff_opt(ii);
    %     single_path = strcat( 'workdir/', Study_paths{ii,1}, '/', Study_paths{ii,2}, '/opt/');
    %     load ( strcat(single_path, 'VOI.mat'));
    %     mu_eff_iter = mu_eff_opt; % Make a copy of both the mu_eff values and the paths
    %     mu_eff_iter ( ii ) = 0; % Set one of the mu_eff values to 0
    %     mu_eff_iter ( mu_eff_iter == 0 ) = []; % Remove the 0
    %     params_iter.cv.mu_eff_healthy = num2str( mean ( mu_eff_iter )); % Average the training datasets' mu_eff; also make it a string coz the thermal code needs that format.
    %     params_iter.patientID = Study_paths{ii,1}; % Write the patient path information into the params_iter structure
    %     params_iter.UID = Study_paths{ii,2};
    %     params_iter.voi(1:2) = VOI.x;
    %     params_iter.voi(1) = 80;
    %     params_iter.voi(3:4) = VOI.y;
    
    % This section runs the thermal code
%     mu_eff(1) = 0.008;
%     mu_eff(2:20001) = linspace(0.5,10000,20000);
    %mu_eff = linspace(1000,5000,8001);
    total = zeros(length(mu_eff_opt),6);
    for jj = 1:length(mu_eff_opt)
        params_iter.cv.mu_eff_healthy = mu_eff_opt (jj);
        
        [metric, thermal_model, MRTI_crop] = fast_temperature_obj_fxn33 ( params_iter );
        % Column 2 of 'total' is based on conservation of energy (only cares
        % about summation of temperatures in the FOV)
        base_level=ones(size(thermal_model,1),size(thermal_model,2))*37;
        total(jj,1) = mu_eff_opt(jj);
        total(jj,2) = sum(sum(thermal_model))-sum(sum(base_level));  
        total(jj,3) = sum(sum(thermal_model))-sum(sum(MRTI_crop)); % No need to subtract base_level (it's in both quantities already, so the difference cancels it)
        total(jj,4) = max(max(thermal_model));
        total(jj,5) = metric;
        %Comment this out (until 'end') if running a long sensitivity
        %study, or maybe just exclude the figure displays.
        model_deg57 = zeros( size(thermal_model,1), size(thermal_model,2) );
        MRTI_deg57 = zeros( size(MRTI_crop,1), size(MRTI_crop,2) );
        model_deg57 = thermal_model >= 57;
        MRTI_deg57 = MRTI_crop >= 57;
        n_model = sum(sum( model_deg57 ));
        n_MRTI = sum(sum( MRTI_deg57 ));
        intersection = model_deg57 + MRTI_deg57;
        intersection = intersection > 1;
        n_intersection = sum(sum( intersection ));
        dice_values = 2*n_intersection / (n_model + n_MRTI) ;
        total(jj,6) = dice_values;
%         figure(1); imagesc(MRTI_crop,[30 80]);
%         figure(2);imagesc(thermal_model, [30 80]);
%         figure(3);imagesc(intersection);
    end
end

% mean_dice = mean ( dice_values ); % Average Dice values
% std_dice = std ( dice_values ); % Stardard deviation of Dice values.
% test_statistic = ( mean_dice - 0.7 ) ./ ( std_dice ./ sqrt ( n_patients )); % Calculate the test_statistic
% p_value = tcdf ( test_statistic , n_patients , 'upwards' );
% H0.p_value = p_value;
% H0.result = 1;
% if p_value <= 0.05
%     H0.result = 0 ;
% else
%     H0.result = 1;
% end

end