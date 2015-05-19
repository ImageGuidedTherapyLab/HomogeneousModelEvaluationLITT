function [FOV_pass] = Patent_LOOCV_figures ( Study_paths, mu_eff_opt, alpha_opt, best_iter, opttype, Matlab_flag, pass_indices, HiLoFlag );

% Make the LOOCV iteration system

n_pass = length(pass_indices);
% n_patients = 1;
%path_iter = cell(1,2);
FOV_pass=zeros(n_pass,2);
%jj=1;
for ii = 1:n_pass
    cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
    disp( fprintf('%d of %d\n',ii,n_pass) );
    
    % Set up LOOCV for mu_eff
    mu_eff_iter = mu_eff_opt; % Make a copy of both the mu_eff values and the paths
    mu_eff_iter ( pass_indices(ii) ) = []; % Remove the 0
    mu_eff_iter = mean ( mu_eff_iter );           % Add alpha to LOOCV here; write *.in.* file; run brainsearch.py; brainsearch.py makes *.out.* file; read in *.out.* file
    
    path_base = strcat( 'workdir/',Study_paths{pass_indices(ii),1},'/',Study_paths{pass_indices(ii),2},'/opt/optpp_pds.',opttype,'.in.1.mat');
    aaa = load(path_base);
    aaa.inputdatavars.cv.mu_eff_healthy = num2str( mu_eff_iter);
    [~,~, tmap,MRTI] =  temperature_obj_fxn ( aaa.inputdatavars, 10 );
    
    aaa.inputdatavars.voi=double(aaa.inputdatavars.voi);
    FOV_pass(ii,1) = aaa.inputdatavars.spacing(1).*(aaa.inputdatavars.voi(2)-aaa.inputdatavars.voi(1));
    FOV_pass(ii,2) = aaa.inputdatavars.spacing(2).*(aaa.inputdatavars.voi(4)-aaa.inputdatavars.voi(3));
    %jj=jj+1;
    model_deg_threshold = tmap >= 57;
    MRTI_deg_threshold = MRTI >= 57;
    n_model = sum(sum( model_deg_threshold ));
    n_MRTI = sum(sum( MRTI_deg_threshold ));
    intersection = model_deg_threshold + 2*MRTI_deg_threshold;

    cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests/Patent_figs
    
    if HiLoFlag == 0
    handle1=figure; imagesc(tmap, [30 80]); colorbar;
    print('-dpng','-r200', strcat('dmg3SS',aaa.inputdatavars.UID,'.png'));
    handle2=figure; imagesc(MRTI, [30 80]); colorbar;
    print('-dpng','-r200', strcat('dmg3MRTI',aaa.inputdatavars.UID,'.png'));
    handle3=figure; imagesc(intersection,[0,3]); colorbar;
    print('-dpng','-r200',strcat('dmg3Int',aaa.inputdatavars.UID,'.png'));
    
    else 
    handle1=figure; imagesc(tmap, [30 80]); colorbar;
    print('-dpng','-r200', strcat('dmg7SS',aaa.inputdatavars.UID,'.png'));
    handle2=figure; imagesc(MRTI, [30 80]); colorbar;
    print('-dpng','-r200', strcat('dmg7MRTI',aaa.inputdatavars.UID,'.png'));
    handle3=figure; imagesc(intersection,[0,3]); colorbar;
    print('-dpng','-r200',strcat('dmg7Int',aaa.inputdatavars.UID,'.png'));
    end
    clear aaa

%     params_iter = load( 'TmpDataInput.mat' ); % Read in one dakota.in file to find the constant parameters
%     params_iter.cv.mu_a = mu_a_iter;
%     params_iter.cv.mu_eff_healthy = num2str( mu_eff_iter ); % Average the training datasets' mu_eff; also make it a string coz the thermal code needs that format.   
%     % This section runs the thermal code
%     [metric, dice_iter, thermal_model, MRTI_crop] = fast_temperature_obj_fxn_sanity ( params_iter, 1 );
%     dice_values (ii) = dice_iter ;
%     clear mu_eff_iter;
%     disp(python_command );
%     %evalc(python_command);
% 
%     % the  result_file  is written by the python command
%     metrics = dlmread(result_file );
%     l2diff          = metrics (1)
%     dice_values(ii) = metrics (2)
    
    % TODO move to your kernel
    %params_iter = load( 'TmpDataInput.mat' ); % Read in one dakota.in file to find the constant parameters
    %%single_path = strcat( 'workdir/', Study_paths{ii,1}, '/', Study_paths{ii,2}, '/opt/');
    %mu_s_p = params_iter.cv.mu_s * ( 1 - params_iter.cv.anfact );
    %mu_a_iter = (-3*mu_s_p + sqrt( 9*mu_s_p^2 + 12 * mu_eff_iter^2))/6; % Used quadratic equation to solve for mu_a in terms of g, mu_s, and mu_eff
    %params_iter.cv.mu_a = mu_a_iter;
    %params_iter.cv.mu_eff_healthy = num2str( mu_eff_iter ); % Average the training datasets' mu_eff; also make it a string coz the thermal code needs that format.
   
    %% This section runs the thermal code
    %[metric, thermal_model, MRTI_crop] = fast_temperature_obj_fxn_sanity ( params_iter, 1 );
    %model_deg57 = zeros( size(thermal_model,1), size(thermal_model,2) );
    %MRTI_deg57 = zeros( size(MRTI_crop,1), size(MRTI_crop,2) );
    %model_deg57 = thermal_model >= 57;
    %MRTI_deg57 = MRTI_crop >= 57;
    %n_model = sum(sum( model_deg57 ));
    %n_MRTI = sum(sum( MRTI_deg57 ));
    %intersection = model_deg57 + MRTI_deg57;
    %intersection = intersection > 1;
    %n_intersection = sum(sum( intersection ));
    %dice_values (ii) = 2*n_intersection / (n_model + n_MRTI) ;
    %clear mu_eff_iter;

end


end
