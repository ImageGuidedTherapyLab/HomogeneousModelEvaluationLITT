% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [total, dice, hd] = temperature_obj_fxn_GPU_perf ( inputdatavars, sources, mu_eff_list, w_perf );
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
mu_length = length(mu_eff_list);
w_length = length(w_perf);
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
%g_anisotropy = str2num(inputdatavars.cv.gamma_healthy);
g_anisotropy = inputdatavars.cv.anfact;
%mu_a = inputdatavars.cv.mu_a;
% mu_s = inputdatavars.cv.mu_s;
% mu_eff = str2num(inputdatavars.cv.mu_eff_healthy); % This should typically be str2num
%k_cond = inputdatavars.cv.k_0;
alpha = str2num(inputdatavars.cv.alpha_healthy);
rho = inputdatavars.cv.rho;
c_p = str2num(inputdatavars.cv.c_p_healthy);
c_blood = str2num(inputdatavars.cv.c_blood_healthy);
%w_perf = inputdatavars.cv.w_0;
geometry.x_disp = str2num(inputdatavars.cv.x_displace);
geometry.y_disp = str2num(inputdatavars.cv.y_displace);
geometry.z_disp = str2num(inputdatavars.cv.z_displace);
geometry.x_rot  = str2num(inputdatavars.cv.x_rotate);
geometry.y_rot  = str2num(inputdatavars.cv.y_rotate);
geometry.z_rot  = str2num(inputdatavars.cv.z_rotate);

% geometry.x_disp = 0;
% geometry.y_disp = 0;
% geometry.z_disp = 0;

if inputdatavars.fileID == 1
    geometry.x_disp = 0;
    geometry.y_disp = 0;
    geometry.z_disp = 0;
%     geometry.x_rot  = 0;
%     geometry.y_rot  = 0;
%     geometry.z_rot  = 0;
else
    aaa = load( strcat( patient_opt_path, 'optpp_pds.', inputdatavars.opttype, '.in.1.mat') );
    
    geometry.x_disp = geometry.x_disp - str2num(aaa.inputdatavars.cv.x_displace);
    geometry.y_disp = geometry.y_disp - str2num(aaa.inputdatavars.cv.y_displace);
    geometry.z_disp = geometry.z_disp - str2num(aaa.inputdatavars.cv.z_displace);
%     geometry.x_rot  = geometry.x_rot - str2num(aaa.inputdatavars.cv.x_rotate);
%     geometry.y_rot  = geometry.y_rot - str2num(aaa.inputdatavars.cv.y_rotate);
%     geometry.z_rot  = geometry.z_rot - str2num(aaa.inputdatavars.cv.z_rotate);
end

% mu_s_p = mu_s * ( 1 - g_anisotropy );
% mu_a = (-3*mu_s_p + sqrt( 9*mu_s_p^2 + 12 * mu_eff^2))/6;

k_cond = alpha * rho * c_p;

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
%mod_point.z = 1;
mod_point.z = 1;
mod_point.z_subslice = 5;

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
% scaling.x = 1;
% scaling.y = 1;
scaling.x = 5;
scaling.y = scaling.x;
scaling.z = 1;

% Build the domain
[dom,~,~]=modeled_domain(FOV,matrix,scaling,mod_point);

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
    
% source.laser = zeros(source.n,3);
% 
% for jj = 1:source.n
%     
%     source.laser(jj,:) = [ source.base(jj)+geometry.x_disp, geometry.y_disp, geometry.z_disp ]; % Columns are for x,y,z displacement; Rows are for different sources
%     
% end
% 
% source.laser = source.laser*rodrigues( [geometry.x_rot*pi/180,(geometry.y_rot+90)*pi/180,geometry.z_rot*pi/180]);

source.laser = zeros(source.n,2);
theta = (geometry.z_rot)*pi/180;
rot_matrix2D = [ cos(theta), -sin(theta); sin(theta), cos(theta) ];
for jj = 1:source.n
    
    source.laser(jj,:) = [ source.base(jj)+geometry.x_disp, geometry.y_disp ]; % Columns are for x,y,z displacement; Rows are for different sources
    source.laser(jj,:) = rot_matrix2D * source.laser(jj,:)';
end
z_dim = zeros(source.n,1)+geometry.z_disp;
source.laser = [source.laser z_dim];

% Run the Bioheat model with the unique powers, and then scale it to MRTI
spacing.x = spacing.x/scaling.x;
spacing.y = spacing.y/scaling.y;
spacing.z = spacing.z/(scaling.z * dom.z_subslice);



[tmap_unique] = Human_GPU_perf ( power_log,spacing,scaling,mod_point,source,w_perf,k_cond,g_anisotropy,mu_eff_list,probe_u,robin_co,c_blood);


%[tmap_unique]=Bioheat1D_SSS( power_log,dom,source,w_perf,k_cond,g_anisotropy,mu_a,mu_s,probe_u,robin_co,c_blood);

tmap_unique=tmap_unique+37;
tmap_model_scaled_to_MRTI = imresize (tmap_unique , 1/scaling.x); % Set the model's spacing to MRTI
 

% Change directory and load the temperature from VTK
cd (patient_MRTI_path);
%MRTI_crop = readVTK_one_time22('roimrtidose.heating', inputdatavars.maxheatid);   % This 'vtkNumber' should b
MRTI = readVTK_one_time('temperature', inputdatavars.maxheatid);   % This 'vtkNumber' should b
% Go to the workdir
cd (patient_opt_path);

% Crop the MRTI using the VOI.
% This is the VOI.x and VOI.y swapped for ParaView
MRTI ( 1:(VOI.y(1)-1), :, : ) = 0;  % The -1 and +1 make sure the VOI indices aren't cut
MRTI ( (VOI.y(2)+1):end, :, : )  = 0;
MRTI ( :, 1:(VOI.x(1)-1), : ) = 0;
MRTI ( :,(VOI.x(2)+1):end, : ) = 0;


MRTI_crop = MRTI( (VOI.y(1) ):(VOI.y(2) ) , (VOI.x(1) ):(VOI.x(2) ) ); % Set the cropped region
MRTI_crop = permute( MRTI_crop, [2 1 3]);
%MRTI_crop22 = MRTI_hottest_DF2( (VOI.y(1) ):(VOI.y(2) ) , (VOI.x(1) ):(VOI.x(2) ) );
% MRTI_crop = MRTI_hottest( (VOI.y(1)-1):(VOI.y(2)+1) , (VOI.x(1)-1):(VOI.x(2)+1) ); % Set the cropped region with odd border

% This is the MATLAB native
%MRTI_crop = MRTI_hottest( (VOI.x(1)-1):(VOI.x(2)+1) , (VOI.y(1)-1):(VOI.y(2)+1) ); % Set the cropped region

% Make the metric

%%%%% There is an error here because ParaView and MATLAB don't agree on
%%%%% row- vs column- major schemes. My solution is to switch VOI.x with
%%%%% VOI.y
cd (path22);
%temperature_diff = tmap_model_scaled_to_MRTI - MRTI_crop ( 2:(end-1), 2:(end-1));
%temperature_diff = tmap_model_scaled_to_MRTI - MRTI_crop;

%temperature_diff = tmap_model_scaled_to_MRTI - MRTI_big;
total = zeros(mu_length,6);
total(:,1) = mu_eff_list;
% mu_eff   L2     weird metric     model-base   max(model)

% isotherms = 51:65
dice = zeros(mu_length,15);
hd = dice;


base_level= sum( sum( ones(size(MRTI_crop,1),size(MRTI_crop,2))*37 ));
for ii = 1:mu_length
    
    % Dice
    for kk=1:15
        model_deg_threshold = tmap_model_scaled_to_MRTI(:,:,ii) >= ( 50 + kk);
        [mod_row, mod_column] = find( model_deg_threshold ==1);
        mod_list = [(mod_row .* inputdatavars.spacing(1)) (mod_column .* inputdatavars.spacing(2))];
        MRTI_deg_threshold = MRTI_crop >= ( 50 + kk);
        [MRTI_row, MRTI_column] = find( MRTI_deg_threshold ==1);
        MRTI_list = [(MRTI_row .* inputdatavars.spacing(1)) (MRTI_column .* inputdatavars.spacing(2))];
        n_model = sum( sum( model_deg_threshold ));
        n_MRTI =  sum( sum( MRTI_deg_threshold  ));
        intersection = model_deg_threshold + MRTI_deg_threshold;
        intersection = intersection > 1;
        n_intersection = sum( sum( intersection ));
        dice(ii,kk) = 2 * n_intersection / ( n_model + n_MRTI );
        if isempty(mod_list)==1
            hd(ii,kk) = 1;
        else
            hd(ii,kk) = HausdorffDist( mod_list, MRTI_list);
        end
    end
    
    % L2
    temperature_diff = tmap_model_scaled_to_MRTI(:,:,ii) - MRTI_crop;
    total(ii,2) = norm ( temperature_diff , 1 ) ; % L1 norm
    
    total(ii,3) = norm ( temperature_diff , 2) ; % L2 norm
    
    total(ii,4) = norm ( temperature_diff , inf); %L_inf norm
    
    % model-base
    total(ii,5) = sum(sum(sum( tmap_model_scaled_to_MRTI(:,:,ii) ))) - base_level;
    
    % max temp
    total(ii,6) = max(max(max( tmap_model_scaled_to_MRTI(:,:,ii) )));


end

end