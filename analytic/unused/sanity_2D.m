% This is a sanity script for the 2-D case
clear
close all
cool_case = 1;
ind_case = 2;
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

% 1D, mu_eff = 180; cooling u0 = 21 C
opt_type = 'sanity_1D' ;
Study_paths {1,1} = 'Study0035';
Study_paths {1,2} = '0530';
if ind_case ==1;
    Study_paths {1,3} = '1';
elseif ind_case ==2;
    Study_paths {1,3} = '43';
end
param_file = strcat( 'workdir/', Study_paths{1,1}, '/', Study_paths{1,2}, '/opt/', opt_type, '.in.', Study_paths{1,3} );
python_command = strcat( 'unix(''python ./brainsearch.py --param_file ./', param_file, ''')');   % unix(''python test_saveFile.py'')
evalc(python_command);

aa = load( 'TmpDataInput.mat' );
aa.cv.body_temp = str2num( aa.cv.body_temp );
Pwr = 12;
R1 = 0.0007;
R2 = 1;
if cool_case == 1;

elseif cool_case ==2;
    aa.cv.probe_init = '37';
end

% Abbrev Run
model = zeros(36,31,16);

% for ii = 1:16
%     [~, model(:,:,ii),MRTI_crop] = fast_temperature_obj_fxn_sanity (aa,ii);
% end

[~, model,MRTI_crop] = fast_temperature_obj_fxn_sanity (aa,3);

model_deg57 = model >= 57;
MRTI_deg57 = MRTI_crop >= 57;
n_model = sum(sum( model_deg57));
n_MRTI = sum(sum( MRTI_deg57 ));
intersection = model_deg57(:,:,1) + MRTI_deg57;
intersection = intersection > 1;
n_intersection = sum(sum( intersection ));
dice_values = 2*n_intersection / (n_model + n_MRTI) ;

% Normal Run
% [~,model1,~]  = fast_temperature_obj_fxn_sanity ( aa,1 );
% [~,model3,~]  = fast_temperature_obj_fxn_sanity ( aa,3 );
% [~,model5,~]  = fast_temperature_obj_fxn_sanity ( aa,5 );
% [~,model11,~] = fast_temperature_obj_fxn_sanity ( aa,11 );
% [~,model14,~] = fast_temperature_obj_fxn_sanity ( aa,14 );
% [~,model15,~] = fast_temperature_obj_fxn_sanity ( aa,15 );
% [~,model16,~] = fast_temperature_obj_fxn_sanity ( aa,16);
% [~,model27,~] = fast_temperature_obj_fxn_sanity ( aa,27 );
% 
% figure;imagesc(model1);
% figure;imagesc(model3);
% figure;imagesc(model5);
% figure;imagesc(model11);
% figure;imagesc(model27);
% 
% hot_spot = zeros(5,1);
% hot_spot(1) = model1(18,16);
% hot_spot(2) = model3(18,16);
% hot_spot(3) = model5(18,16);
% hot_spot(4) = model11(18,16);
% hot_spot(5) = model27(18,16);
% hot_spot