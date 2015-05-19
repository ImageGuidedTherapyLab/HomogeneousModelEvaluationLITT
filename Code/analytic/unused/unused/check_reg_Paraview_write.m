% This function is meant to find the registration for the ROIs. The first
% argin is reg and should be an array [ x y z ]. But note that the VOI.x,y indexices are swapped. The pwd should be:
% /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/

% function [ center, VOI ] = check_reg_write (  reg,time,crop_big,crop_small );
function [ center, VOI ] = check_reg_Paraview_write (  VOI_pre, path22, work_dir, pathpt );

% Get to the laser directory
total_laser_path = strcat ( path22, pathpt); % This path needs mofication by strrep
total_laser_path = strrep ( total_laser_path, '/workdir', '');
total_laser_path = strrep ( total_laser_path, 'opt', 'laser');

cd (total_laser_path);


% Import the VTK header info
FOV_import = csvimport ( 'FOV.csv' );
fov = FOV_import ( 2, : );
fov = cell2mat (fov);

% Define the domain and scaling
mod_point.x=171;  % x position on image
mod_point.y=171;
mod_point.z=1;

%mod_center=[floor(mod_point.y/2) floor(mod_point.x/2) 1]; % [y x z]
% Record the vtk header here
matrix.x = fov(1);
matrix.y = fov(2);
matrix.z = fov(3);

spacing.x = fov(4);
spacing.y = fov(5);
spacing.z = fov(6);

FOV.x = matrix.x * spacing.x;
FOV.y = matrix.y * spacing.y;
FOV.z = matrix.z * spacing.z;

% For now, x and y scaling must be equal; z = 1
scaling.x=1;
scaling.y=1;
scaling.z=1;

% Build the domain
[~,MRTI_pix,~]=modeled_domain(FOV,matrix,scaling,mod_point);

reg = VOI_pre.center;

reg(1) = reg(1) * MRTI_pix.x; % Change registration and divide by MRTI pixel size.
reg(2) = reg(2) * MRTI_pix.y;
reg(3) = reg(3) * MRTI_pix.z;
    
%Find the center point
center.x = round ( reg (1) ./ MRTI_pix.x );
center.y = round ( reg (2) ./ MRTI_pix.y );
center.z = round ( reg (3) ./ MRTI_pix.z );

center_in_meters = reg (1:3);

% Display a volume of interest, VOI

% The commented stuff is if the input for registration is [ x y z ]
% VOI_small.x = [(round ( (reg (1) ./ MRTI_pix.x) - crop_small )) (round ( (reg (1) ./ MRTI_pix.x) + crop_small ))];
% VOI_small.y = [(round ( (reg (2) ./ MRTI_pix.y) - crop_small )) (round ( (reg (2) ./ MRTI_pix.y) + crop_small ))];
% VOI_small.z = [(round ( (reg (3) ./ MRTI_pix.z) - crop_small )) (round ( (reg (3) ./ MRTI_pix.z) + crop_small ))];

% VOI.x = [(round ( (reg (1) ./ MRTI_pix.x) - crop_big )) (round ( (reg (1) ./ MRTI_pix.x) + crop_big ))];
% VOI.y = [(round ( (reg (2) ./ MRTI_pix.y) - crop_big )) (round ( (reg (2) ./ MRTI_pix.y) + crop_big ))];
% VOI.z = [(round ( (reg (3) ./ MRTI_pix.z) - crop_big )) (round ( (reg (3) ./ MRTI_pix.z) + crop_big ))];

% VOI.x = [(round ( (reg (1) ./ MRTI_pix.x) - crop_big )) (round ( (reg (1) ./ MRTI_pix.x) + crop_big ))];
% VOI.y = [(round ( (reg (2) ./ MRTI_pix.y) - crop_big )) (round ( (reg (2) ./ MRTI_pix.y) + crop_big ))];
% VOI.z = [(round ( (reg (3) ./ MRTI_pix.z) - crop_big )) (round ( (reg (3) ./ MRTI_pix.z) + crop_big ))];

VOI.x = VOI_pre.x;
VOI.y = VOI_pre.y;
VOI.z = VOI_pre.x;

% VOI_small.x = [(round ( (reg (1) ./ MRTI_pix.x) - crop_small )) (round ( (reg (1) ./ MRTI_pix.x) + crop_small ))];
% VOI_small.y = [(round ( (reg (2) ./ MRTI_pix.y) - crop_small )) (round ( (reg (2) ./ MRTI_pix.y) + crop_small ))];
% VOI_small.z = [(round ( (reg (3) ./ MRTI_pix.z) - crop_small )) (round ( (reg (3) ./ MRTI_pix.z) + crop_small ))];

% cd ../matlab
% load 'referenceBased.mat'
% MRTI = (temperature);

cd ../vtk/referenceBased
MRTI = readVTK_SJF('temperature',400);

% Crop the displayed VOI
MRTI_crop = MRTI;                            % Start with original MRTI
MRTI_crop ( 1 : (VOI.x(1)-1)   , :,:,:) = 0; % Crop top
MRTI_crop ( (VOI.x(2)+1) : end , :,:,:) = 0; % Overwrite MRTI_crop, crop bottom
MRTI_crop ( : ,  1 : (VOI.y(1)-1) ,:,:) = 0; % Crop left
MRTI_crop ( : , (VOI.y(2)+1) : end,:,:) = 0; % Crop right

% MRTI_crop_small = MRTI;                                  % Start with original MRTI
% MRTI_crop_small ( 1 : (VOI_small.x(1)-1)   , :,:,:) = 0; % Crop top
% MRTI_crop_small ( (VOI_small.x(2)+1) : end , :,:,:) = 0; % Overwrite MRTI_crop, crop bottom
% MRTI_crop_small ( : ,  1 : (VOI_small.y(1)-1) ,:,:) = 0; % Crop left
% MRTI_crop_small ( : , (VOI_small.y(2)+1) : end,:,:) = 0; % Crop right

% figure(1); imagesc ( MRTI (:,:,VOI_pre.time) , [40 80] );
% figure(2); imagesc ( MRTI_crop (:,:,VOI_pre.time) , [40 80] );
% figure(3); imagesc ( MRTI_crop_small (:,:,VOI_pre.time) , [40 80] );

VOI.center_in_meters = reg;
VOI.center_in_pix = [ center.x center.y center.z ];
VOI.time = VOI_pre.time;
center;
center_in_meters;
VOI;

total_opt_path = strcat ( work_dir, pathpt); % This path needs mofication by strrep
cd ( total_opt_path );
save ('VOI.mat' , 'VOI') ;
end