% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [dice] = temperature_GPU_LOOCV_rerun ( inputdatavars, sources, max_phys_sz, mu_eff_list, w_perf, k_cond, choice );
%[total, dice, hd, mutual_threshold, false_pix] = temperature_obj_fxn_GPU_choice ( inputdatavars, sources, mu_eff_list, w_perf, k_cond, choice );
% Record the working directory
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% Make the path to the patient directory
patientID = strcat ( inputdatavars.patientID, '/', inputdatavars.UID, '/');
patient_opt_path = strcat ( path22, '/workdir/', patientID, 'opt/' ); % The '/workdir/' needs the first backslash coz it uses absolute path

% Make the path to the VTK
patient_MRTI_path = strcat ( 'StudyDatabase/', patientID, 'vtk/referenceBased/' );
%patient_MRTI_path = strcat ( '/tmp/outputs/dakota/',inputdatavars.UID,'/');
% % Read in the power and identify the max power moment.
% power_log = str2num(inputdatavars.powerhistory); % Gets just the numbers. (Al a "Just the facts, ma'am.")


% if choice ==1
%     n_length = length(mu_eff_list);
%     
% elseif choice ==2
%     n_length = length(w_perf);
%     
% elseif choice ==3
%     n_length = length(k_cond);
%     
% elseif choice ==4
%     w_length = length(w_perf);
%     m_length = length(mu_eff_list);
%     
% end
scaling.x = 5;
scaling.y = 5;
scaling.z = 1;
mod_point.z_subslice = 5;


spacing.x = max_phys_sz(1,3)/scaling.x;
spacing.y = max_phys_sz(2,3)/scaling.y;
spacing.z = 0.003;
spacing.z = spacing.z/( scaling.z * mod_point.z_subslice -1 );  % The "-1" is to get the number of intervals correct. Imagine 2 slices separated by .03. There is 1 interval of 0.03 distance, not two intervals of 0.015.

mod_point.x = ceil((max_phys_sz(1,2).* max_phys_sz(1,4) ) / max_phys_sz(1,3) );
mod_point.y = ceil((max_phys_sz(2,2).* max_phys_sz(2,4) ) / max_phys_sz(2,3) );
mod_point.z = 1;

sim_dim.spacing = spacing;
sim_dim.mod_point = mod_point;

g_anisotropy = inputdatavars.cv.anfact;
probe_u = str2num(inputdatavars.cv.probe_init);
robin_co=0; %dummy var
c_blood = str2num(inputdatavars.cv.c_blood_healthy);

% Define the source; source.n must be greater than 1 to make sense. Odd is
% better than even
source.n=sources;
source.length=0.01;  %~0.033 is when n=5 is visible
if source.n > 1;
    source.base = linspace((-source.length/2),(source.length/2),source.n);
    % source.laser = [geometry.x_disp, geometry.y_disp, geometry.z_disp];
else
    source.base = 0;
end

source.laser = zeros(source.n,3);
source.laser(:,1) = source.base;

% Cold bases ('normal')
% power_cold = 0;
% k_cold = 0.527;
% w_cold = 6;
% mu_cold = 200;
% if choice == 1
%     
%     [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%     
% elseif choice == 2
%     
%     w_cold = w_perf;
%     [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%     
% elseif choice == 3
%     
%     k_cold = k_cond;
%     [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%     
% elseif choice == 5
%     choice_cold=2;
%     w_cold = w_perf;
%     [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice_cold);
%     
% end

pwr_hist = str2num(inputdatavars.powerhistory);
for ii = 1:(length(pwr_hist) - 1)  % Write the for loop to iterate through all but the last index of 'pwr_hsitory'
    diff = pwr_hist(ii) - pwr_hist( ii + 1); % Record the difference between neighboring array elements
    if diff >= 0 % If 'diff' is non-negative, it means the times have stopped and the powers are starting.
        break;  % Stop the for-loop.      
    end
end
% Based on immediately previous for-loop, parse the times from powers.
power_log = pwr_hist ( (ii+1):end);
power_log = max(power_log); % Find the maximum power value

clear diff

% Make the VOI; Note the 'inputdatavars.voi' is from ParaView.

% Make the VOI; Note the 'inputdatavars.voi' is from ParaView.
inputdatavars.voi = double( inputdatavars.voi );
VOI.x = inputdatavars.voi(3:4); % The weird index assignment is coz it's from ParaView.
VOI.y = inputdatavars.voi(1:2);
VOI.z = inputdatavars.voi(5:6);

if sum( inputdatavars.UID == '0496' ) ==4
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0402' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 3;
elseif sum( inputdatavars.UID == '0389' ) ==4
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 2;
elseif sum( inputdatavars.UID == '0385' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y - 1;
    inputdatavars.maxheatid = 109;
elseif sum( inputdatavars.UID == '0476' ) ==4
    VOI.x = VOI.x + 5;
    VOI.y = VOI.y - 5;
elseif sum( inputdatavars.UID == '0477' ) ==4
    VOI.x = VOI.x + 3;
    VOI.y = VOI.y - 2;
elseif sum( inputdatavars.UID == '0438' ) ==4
    VOI.x = VOI.x - 1;
    VOI.y = VOI.y + 0;
elseif sum( inputdatavars.UID == '0435' ) ==4
    VOI.x = VOI.x - 3;
    VOI.y = VOI.y + 3;
elseif sum( inputdatavars.UID == '0436' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 6;
    power_log = 10.05;
    inputdatavars.maxheatid = 39;
elseif sum( inputdatavars.UID == '0466' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
    power_log = 12;
elseif sum( inputdatavars.UID == '0468' ) ==4
    VOI.x = VOI.x - 3;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0471' ) ==4
    VOI.x = VOI.x + 3;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0447' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0453' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0451' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0418' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 0;
    %inputdatavars.maxheatid = 85;
elseif sum( inputdatavars.UID == '0409' ) ==4
    VOI.x = VOI.x - 2;
    VOI.y = VOI.y + 11;
elseif sum( inputdatavars.UID == '0414' ) ==4
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0415' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 1;
end
[model_temp] = Human_GPU_choice_sym ( power_log,spacing,scaling,mod_point,source,w_perf,k_cond,g_anisotropy,mu_eff_list,probe_u,robin_co,c_blood,choice);

model_temp = imrotate(model_temp, str2num( inputdatavars.cv.z_rotate),'crop' );
model_temp = model_temp + 37; % Get to body temp 

% Set up MRTI meshgrid
FOV_l(1) = inputdatavars.spacing(1) .* (inputdatavars.voi(2)-inputdatavars.voi(1)+1);
FOV_l(2) = inputdatavars.spacing(2) .* (inputdatavars.voi(4)-inputdatavars.voi(3)+1);
%[MRXq, MRYq] = meshgrid ( inputdatavars.spacing(1): inputdatavars.spacing(1) : FOV_l(1) , inputdatavars.spacing(2): inputdatavars.spacing(2) : FOV_l(2) );
[MRXq, MRYq] = meshgrid ( -(FOV_l(1)-inputdatavars.spacing(1))/2: inputdatavars.spacing(1) : (FOV_l(1)-inputdatavars.spacing(1))/2 , -(FOV_l(2)-inputdatavars.spacing(2))/2: inputdatavars.spacing(2) : (FOV_l(2)-inputdatavars.spacing(2))/2 );

% Set up model meshgrid

sim_dim.spacing = spacing;
sim_dim.mod_point = mod_point;

modFOV_l(1) = sim_dim.spacing.x .*sim_dim.mod_point.x .*sim_dim.mod_point.z_subslice;
modFOV_l(2) = sim_dim.spacing.y .*sim_dim.mod_point.y .*sim_dim.mod_point.z_subslice;
[modX, modY] = meshgrid ( -(modFOV_l(1)-sim_dim.spacing.x)/2: sim_dim.spacing.x : (modFOV_l(1)-sim_dim.spacing.x)/2 , -(modFOV_l(2)-sim_dim.spacing.y)/2: sim_dim.spacing.y : (modFOV_l(2)-sim_dim.spacing.y)/2 );

% Do interpolation
%n_model = size( model_temp,3);
%model_crop = zeros( (inputdatavars.voi(4)-inputdatavars.voi(3)+1), (inputdatavars.voi(2)-inputdatavars.voi(1)+1));
model_crop = qinterp2( modX, modY, model_temp, MRXq, MRYq);
model_crop=model_crop';

cd (path22);
% Change directory and load the temperature from VTK
cd (patient_MRTI_path);


% Make the VOI; Note the 'inputdatavars.voi' is from ParaView.
MRTI = readVTK_one_time('temperature', inputdatavars.maxheatid);   % This 'vtkNumber' should b
cd (patient_opt_path);

% Crop the MRTI using the VOI.
% This is the VOI.x and VOI.y swapped for ParaView
MRTI ( 1:(VOI.y(1)-1), :, : ) = 0;  % The -1 and +1 make sure the VOI indices aren't cut
MRTI ( (VOI.y(2)+1):end, :, : )  = 0;
MRTI ( :, 1:(VOI.x(1)-1), : ) = 0;
MRTI ( :,(VOI.x(2)+1):end, : ) = 0;

MRTI_crop = MRTI( (VOI.y(1) ):(VOI.y(2) ) , (VOI.x(1) ):(VOI.x(2) ) ); % Set the cropped region
%MRTI_crop = permute( MRTI_crop, [2 1 3]);

cd (path22);

MRTI_isotherm = MRTI_crop >= 57;
MRTI_isotherm = imfill( ExtractNLargestBlobs(MRTI_isotherm,1) , 'holes');
n_MRTI =  sum(  MRTI_isotherm(:) );

model_isotherm = model_crop >= 57;
n_model = sum( model_isotherm(:) );

intersection = model_isotherm + MRTI_isotherm;
intersection = intersection > 1;
n_intersection = sum( intersection(:));
dice = 2 * n_intersection / (n_model + n_MRTI);

label_map = model_isotherm + 2 .* MRTI_isotherm;
figure;imagesc(label_map);
close

end

% % Cold bases ('normal')
% power_cold = 0;
% k_cold = 0.527;
% w_cold = 6;
% mu_cold = 200;
% [no_pwr_fig1] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 0;
% k_cold = 1;
% w_cold = 6;
% mu_cold = 200;
% [no_pwr_fig2] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 0;
% k_cold = 0.527;
% w_cold = 12;
% mu_cold = 200;
% [no_pwr_fig3] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% % 1 W cases
% power_cold = 1;
% k_cold = 0.527;
% w_cold = 6;
% mu_cold = 200;
% [lo_pwr_fig1] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 1;
% k_cold = 1;
% w_cold = 6;
% mu_cold = 200;
% [lo_pwr_fig2] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 1;
% k_cold = 0.527;
% w_cold = 12;
% mu_cold = 200;
% [lo_pwr_fig3] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 1;
% k_cold = 0.527;
% w_cold = 6;
% mu_cold = 500;
% [lo_pwr_fig4] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% % 10 W cases
% power_cold = 10;
% k_cold = 0.527;
% w_cold = 6;
% mu_cold = 200;
% [hi_pwr_fig1] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 10;
% k_cold = 1;
% w_cold = 6;
% mu_cold = 200;
% [hi_pwr_fig2] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 10;
% k_cold = 0.527;
% w_cold = 12;
% mu_cold = 200;
% [hi_pwr_fig3] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
%
% power_cold = 10;
% k_cold = 0.527;
% w_cold = 6;
% mu_cold = 500;
% [hi_pwr_fig4] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);