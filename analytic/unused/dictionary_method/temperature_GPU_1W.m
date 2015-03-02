% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [all_opt_fig, no_pwr_fig,sim_dim] = temperature_GPU_1W ( inputdatavars, sources, max_phys_sz, mu_eff_list, w_perf, k_cond, choice );
% Record the working directory
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% Make the path to the patient directory
patientID = strcat ( inputdatavars.patientID, '/', inputdatavars.UID, '/');
patient_opt_path = strcat ( path22, '/workdir/', patientID, 'opt/' ); % The '/workdir/' needs the first backslash coz it uses absolute path

% Make the path to the VTK
patient_MRTI_path = strcat ( 'StudyDatabase/', patientID, 'vtk/referenceBased/' );
%patient_MRTI_path = strcat ( '/tmp/outputs/dakota/',inputdatavars.UID,'/');
% Read in the power and identify the max power moment.
pwr_hist = str2num(inputdatavars.powerhistory); % Gets just the numbers. (Al a "Just the facts, ma'am.")


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
power_cold = 0;
k_cold = 0.527;
w_cold = 6;
mu_cold = 200;
if choice == 1
    
    [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
    
elseif choice == 2
    
    w_cold = w_perf;
    [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
    
elseif choice == 3
    
    k_cold = k_cond;
    [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice);
    
elseif choice == 5
    choice_cold=2;
    w_cold = w_perf;
    [no_pwr_fig] = Human_GPU_choice_sym ( power_cold,spacing,scaling,mod_point,source,w_cold,k_cold,g_anisotropy,mu_cold,probe_u,robin_co,c_blood,choice_cold);
    
end

power_log = 1;
[all_opt_fig] = Human_GPU_choice_sym ( power_log,spacing,scaling,mod_point,source,w_perf,k_cond,g_anisotropy,mu_eff_list,probe_u,robin_co,c_blood,choice);

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