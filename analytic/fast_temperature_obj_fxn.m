% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [metric,tmap_model_scaled_to_MRTI,MRTI_crop] = fast_temperature_obj_fxn ( inputdatavars );
% Record the working directory
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% Make the path to the patient directory
patientID = strcat ( inputdatavars.patientID, '/', inputdatavars.UID, '/');
patient_opt_path = strcat ( path22, '/workdir/', patientID, 'opt/' ); % The '/workdir/' needs the first backslash coz it uses absolute path

% Make the path to the VTK
patient_MRTI_path = strcat ( 'StudyDatabase/', patientID, 'vtk/referenceBased/' );

% Read in the power and identify the max power moment.
pwr_hist = str2num(inputdatavars.powerhistory); % Gets just the numbers. (Al a "Just the facts, ma'am.")
% This for loop identifies the power history's times versus powers for
% later parsing.
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

% Read in the CVs from inputdatavars . Some need str2num coz they were
% written as strings.
probe_u = str2num(inputdatavars.cv.probe_init);
g_anisotropy = str2num(inputdatavars.cv.gamma_healthy);
mu_a = inputdatavars.cv.mu_a;
mu_s = inputdatavars.cv.mu_s;
mu_eff = str2num(inputdatavars.cv.mu_eff_healthy);
k_cond = inputdatavars.cv.k_0;
w_perf = inputdatavars.cv.w_0;
x_disp = str2num(inputdatavars.cv.x_displace);
y_disp = str2num(inputdatavars.cv.y_displace);
z_disp = str2num(inputdatavars.cv.z_displace);
x_rot  = str2num(inputdatavars.cv.x_rotate);
y_rot  = str2num(inputdatavars.cv.y_rotate);
z_rot  = str2num(inputdatavars.cv.z_rotate);

robin_co=0; %dummy var

% Make the VOI; Note the 'inputdatavars.voi' is from ParaView.
VOI.x = double( inputdatavars.voi(3:4)); % The weird index assignment is coz it's from ParaView.
VOI.y = double( inputdatavars.voi(1:2));
VOI.z = double( inputdatavars.voi(5:6));

% Define the domain and scaling
mod_point.x = abs ( VOI.x(1) - VOI.x(2) ) +1;  % x dimension distance
% The '+ 1' is important becasue there are **51** pixels for
% VOI.x(1):VOI.x(2) but only 50 for VOI.x(1)-VOI.x(2)
mod_point.y = abs ( VOI.y(1) - VOI.y(2) ) +1;  % y dimension
mod_point.z = 1;

% Set the matrix for the MRTI
matrix.x = double(inputdatavars.dimensions(2)); % The indexing is to match ParaView, but all images are 256x256x1 so far, so x,y indexing can be swapped.
matrix.y = double(inputdatavars.dimensions(1));
matrix.z = double(inputdatavars.dimensions(3));

% Set the spacing for the MRTI
spacing.x = double(inputdatavars.spacing(2));
spacing.y = double(inputdatavars.spacing(1));
spacing.z = double(inputdatavars.spacing(3));

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

% Run the Bioheat model with the unique powers, and then scale it to MRTI
[tmap_unique]=Bioheat1Dfast( power_log,dom,source,w_perf,k_cond,g_anisotropy,mu_a,mu_s,probe_u,robin_co);
tmap_unique=tmap_unique+37;
tmap_model_scaled_to_MRTI = imresize (tmap_unique , 1/scaling.x); % Set the model's spacing to MRTI
 

% Change directory and load the temperature from VTK
cd (patient_MRTI_path);
MRTI = readVTK_SJF('temperature', pwr_hist(ii));   % This 'vtkNumber' should b

% Go to the workdir
cd (patient_opt_path);

% Crop the MRTI using the VOI.
% This is the VOI.x and VOI.y swapped for ParaView
MRTI ( 1:(VOI.y(1)-1), :, : ) = 0;  % The -1 and +1 make sure the VOI indices aren't cut
MRTI ( (VOI.y(2)+1):end, :, : )  = 0;
MRTI ( :, 1:(VOI.x(1)-1), : ) = 0;
MRTI ( :,(VOI.x(2)+1):end, : ) = 0;

% This is the MATLAB native
% MRTI ( 1:(VOI.x(1)-1), :, : ) = 0;  % The -1 and +1 make sure the VOI indices aren't cut
% MRTI ( (VOI.x(2)+1):end, :, : )  = 0;
% MRTI ( :, 1:(VOI.y(1)-1), : ) = 0;
% MRTI ( :,(VOI.y(2)+1):end, : ) = 0;

% Find the time with the hottest MRTI and then set a variable to that hot
% timepoint.
timepoint_summation = squeeze( sum( sum ( MRTI , 1 ) , 2 ) ); % Sum the temperature spatially
[~ ,max_temperature_timepoint] = max ( timepoint_summation ); % Find the max summed value for temperature
MRTI_hottest = MRTI(:,:,max_temperature_timepoint); % Set the variable

% Crop the MRTI_image
% This is the VOI.x and VOI.y swapped for ParaView
MRTI_crop = MRTI_hottest( (VOI.y(1) ):(VOI.y(2) ) , (VOI.x(1) ):(VOI.x(2) ) ); % Set the cropped region
% MRTI_crop = MRTI_hottest( (VOI.y(1)-1):(VOI.y(2)+1) , (VOI.x(1)-1):(VOI.x(2)+1) ); % Set the cropped region with odd border

% This is the MATLAB native
%MRTI_crop = MRTI_hottest( (VOI.x(1)-1):(VOI.x(2)+1) , (VOI.y(1)-1):(VOI.y(2)+1) ); % Set the cropped region

% Make the metric

%%%%% There is an error here because ParaView and MATLAB don't agree on
%%%%% row- vs column- major schemes. My solution is to switch VOI.x with
%%%%% VOI.y

%temperature_diff = tmap_model_scaled_to_MRTI - MRTI_crop ( 2:(end-1), 2:(end-1));
temperature_diff = tmap_model_scaled_to_MRTI - MRTI_crop;
metric = ( norm ( temperature_diff , 2 ) )^2;
cd (path22);
% %%%%% The remaining code is exclusively for testing / debugging / checking
% %%%%% the registration.
% MRTI_size = size ( MRTI );
% % Resize the tmap_unique model into the same spacing as the MRTI and find
% % the max heating
% aa = imresize (tmap_unique , 1/scaling.x);  
% 
% % This section writes the model data to the upper left of a matrix
% aa_size = size ( aa );
% size_diff=[(MRTI_size(1)-aa_size(1)) (MRTI_size(2)-aa_size(2))];
% upper_left_mod = zeros((size(aa,1)+size_diff(1)),(size(aa,2)+size_diff(2)));
% upper_left_mod(1:size(aa,1),1:size(aa,2)) = aa; % Write the data to the upper left
% 
% % Register the model data to the MRTI
% matched_mod = zeros (MRTI_size(1), MRTI_size(2));
% matched_mod ( VOI.x(1): VOI.x(2), VOI.y(1): VOI.y(2) ) = upper_left_mod( 1:aa_size(1) , 1:aa_size(2) );  % Write the data to the correct region
% 
% % This is useful for confirming registration, but not for running with
% % DAKMATLAB.
% 
% figure(1); imagesc(tmap_model_scaled_to_MRTI );
% figure(2); imagesc(MRTI_crop , [30 80]);
% figure(3); imagesc(temperature_diff );
% figure(4); imagesc(matched_mod);
% figure(5); imagesc(MRTI(:,:,max_temperature_timepoint), [30 80]);

end