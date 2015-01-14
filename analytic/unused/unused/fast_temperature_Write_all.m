% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [metric] = fast_temperature_Write_all ( path22, pathpt, iteration );
cd( path22);
patient_path = pathpt;
patient_opt_path = strcat( path22, patient_path);
cd( patient_opt_path);
input_param = 'optpp_pds.in.';
index = num2str(iteration);
input_filename = strcat( input_param, index);
input_filename = strcat( input_filename, '.mat' );

load(input_filename);
% aaa=strtrim(aaa);
% aaa=regexp(aaa,'\s+','split'); % Separate the label from the data into two columns.

% Write every string as a number
probe_u = str2num(probe_init);
g_anisotropy = str2num(anfact_healthy);
mu_a = str2num(mu_a_healthy);
mu_s = str2num(mu_s_healthy);
k_cond = str2num(k_0_healthy);
w_perf = str2num(w_0_healthy);
% x_disp = str2num(aaa{8}{1});
% y_disp = str2num(aaa{9}{1});
% z_disp = str2num(aaa{10}{1});
% x_rot = str2num(aaa{11}{1});
% y_rot = str2num(aaa{12}{1});
% z_rot = str2num(aaa{13}{1});

robin_co=0; %dummy var

% Load the recorded power
power_log = load ('time_then_power.csv');

% Change directory to load MATLAB registration data.
dose_path = strcat( '/FUS4/data2/BioTex/BrainNonMDA/processed' , patient_path );
dose_path = strrep ( dose_path, '/workdir' , '' );
dose_path = strrep ( dose_path, '/opt' , '' );
dose_path = strcat( dose_path , '/matlab' );
cd ( dose_path );
load ( 'VOI.mat' );

% Define the domain and scaling
mod_point.x = abs ( VOI.x(1) - VOI.x(2) ) +1;  % x dimension distance
% The '+ 1' is important becasue there are **51** pixels for
% VOI.x(1):VOI.x(2) but only 50 for VOI.x(1)-VOI.x(2)
mod_point.y = abs ( VOI.y(1) - VOI.y(2) ) +1;  % y dimension
mod_point.z = 1;

% Import the VTK header info to determine the MRTI image data dimensions
path_append = strrep ( patient_opt_path, '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/', '');
path_append = strrep ( path_append, 'opt', '');
FOV_path = strcat( '/FUS4/data2/BioTex/BrainNonMDA/processed/' , path_append );
FOV_path = strcat( FOV_path , 'laser' );
cd ( FOV_path );
FOV_import = csvimport ( 'FOV.csv' );
fov = FOV_import ( 2, : );
fov = cell2mat (fov);

% Set the matrix for the MRTI
matrix.x = fov(1);
matrix.y = fov(2);
matrix.z = fov(3);

% Set the spacing for the MRTI
spacing.x = fov(4);
spacing.y = fov(5);
spacing.z = fov(6);

% Set the FOV for the MRTI
FOV.x = matrix.x * spacing.x;
FOV.y = matrix.y * spacing.y;
FOV.z = matrix.z * spacing.z;

% For now, x and y scaling must be equal; z = 1
scaling.x = 1;
scaling.y = 1;
scaling.z = 1;

% Build the domain
[dom,~,~]=modeled_domain(FOV,matrix,scaling,mod_point);

% Define the source; source.n must be greater than 1 to make sense. Odd is
% better than even
source.n=5;
source.length=0.01;  %~0.033 is when n=5 is visible
source.laser=linspace((-source.length/2),(source.length/2),source.n);

% Only use the hottest time point
power_log = max( power_log(:,2));

% Run the Bioheat model with the unique powers, and then scale it to MRTI
[tmap_unique]=Bioheat1Dfast( power_log,dom,source,w_perf,k_cond,g_anisotropy,mu_a,mu_s,probe_u,robin_co);
tmap_unique=tmap_unique+37;
tmap_model_scaled_to_MRTI = imresize (tmap_unique , 1/scaling.x); % Set the model's spacing to MRTI
 

% Change directory and load the temperature from VTK
cd ../vtk/referenceBased
MRTI = readVTK_SJF('temperature',161);

% Crop the MRTI using the VOI.
MRTI ( 1:(VOI.x(1)-1), :, : ) = 0;  % The -1 and +1 make sure the VOI indices aren't cut
MRTI ( (VOI.x(2)+1):end, :, : )  = 0;
MRTI ( :, 1:(VOI.y(1)-1), : ) = 0;
MRTI ( :,(VOI.y(2)+1):end, : ) = 0;

% Find the time with the hottest MRTI and then set a variable to that hot
% timepoint.
timepoint_summation = squeeze( sum( sum ( MRTI , 1 ) , 2 ) ); % Sum the temperature spatially
[~ ,max_temperature_timepoint] = max ( timepoint_summation ); % Find the max summed value for temperature
MRTI_hottest = MRTI(:,:,max_temperature_timepoint); % Set the variable

% Crop the MRTI_image
MRTI_crop = MRTI_hottest( (VOI.x(1)-1):(VOI.x(2)+1) , (VOI.y(1)-1):(VOI.y(2)+1) ); % Set the cropped region

% Make the metric
temperature_diff = tmap_model_scaled_to_MRTI - MRTI_crop ( 2:(end-1), 2:(end-1));
metric = ( norm ( temperature_diff , 2 ) )^2;

% %%%%% The remaining code is exclusively for testing / debugging / checking
% %%%%% the registration.
MRTI_size = size ( MRTI );
% Resize the tmap_unique model into the same spacing as the MRTI and find
% the max heating
aa = imresize (tmap_unique , 1/scaling.x);  

% This section writes the model data to the upper left of a matrix
aa_size = size ( aa );
size_diff=[(MRTI_size(1)-aa_size(1)) (MRTI_size(2)-aa_size(2))];
upper_left_mod = zeros((size(aa,1)+size_diff(1)),(size(aa,2)+size_diff(2)));
upper_left_mod(1:size(aa,1),1:size(aa,2)) = aa; % Write the data to the upper left

% Register the model data to the MRTI
matched_mod = zeros (MRTI_size(1), MRTI_size(2));
matched_mod ( VOI.x(1): VOI.x(2), VOI.y(1): VOI.y(2) ) = upper_left_mod( 1:aa_size(1) , 1:aa_size(2) );  % Write the data to the correct region

% This is useful for confirming registration, but not for running with
% DAKMATLAB.
% HotHot = zeros( size(matched_mod) ); % Must initialize HotHot for the isotherm contour in the next line
% HotHot (matched_mod > 57 ) = 1;
% figure(1); imagesc(tmap_model_scaled_to_MRTI );
% figure(2); imagesc(MRTI_crop , [30 80]);
% figure(3); imagesc(temperature_diff );
% figure(4); imagesc(matched_mod);
% figure(5); imagesc(MRTI(:,:,max_temperature_timepoint), [30 80]);
% figure(6); imagesc( HotHot);

% This section records the VTK files
cd( patient_opt_path);

header.ImagePositionPatient(1) = 0;
header.ImagePositionPatient(2) = 0;
header.ImagePositionPatient(3) = 0;
header.PixelSpacing(1) = spacing.x;
header.PixelSpacing(2) = spacing.y;
header.SliceThickness  = spacing.z;

writeVTK_SJF( matched_mod, 'model_temperature', header );

cd (path22);

end