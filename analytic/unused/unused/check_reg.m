% THis function is meant to find the registration for the ROIs.
% THe pwd should be
% /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/laser

% Once the right center vector is found, make sure they're are entered [ y x
% z ] for check_reg_write.m

function [ center, VOI ] = check_reg ( path22 ,time,crop_big,crop_small );
cd ( path22 );
reg_laser = csvimport ( 'reg_laser.csv' );
reg = cell2mat( reg_laser(2,:) );           % Grabs the numeric portion of the registration cell array
%clear reg;

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

%reg = [(104.*MRTI_pix.x) (128 .*MRTI_pix.y)   (0.*MRTI_pix.z)]; % Change registration and divide by MRTI pixel size.
    
%Find the center point
center.x = round ( reg (1) ./ MRTI_pix.x );
center.y = round ( reg (2) ./ MRTI_pix.y );
center.z = round ( reg (3) ./ MRTI_pix.z );

center_in_meters = reg (1:3);

%Display a volume of interest, VOI
VOI.x = [(round ( (reg (1) ./ MRTI_pix.x) - crop_big )) (round ( (reg (1) ./ MRTI_pix.x) + crop_big ))];
VOI.y = [(round ( (reg (2) ./ MRTI_pix.y) - crop_big )) (round ( (reg (2) ./ MRTI_pix.y) + crop_big ))];
VOI.z = [(round ( (reg (3) ./ MRTI_pix.z) - crop_big )) (round ( (reg (3) ./ MRTI_pix.z) + crop_big ))];

VOI_small.x = [(round ( (reg (1) ./ MRTI_pix.x) - crop_small )) (round ( (reg (1) ./ MRTI_pix.x) + crop_small ))];
VOI_small.y = [(round ( (reg (2) ./ MRTI_pix.y) - crop_small )) (round ( (reg (2) ./ MRTI_pix.y) + crop_small ))];
VOI_small.z = [(round ( (reg (3) ./ MRTI_pix.z) - crop_small )) (round ( (reg (3) ./ MRTI_pix.z) + crop_small ))];

cd ../vtk/referenceBased
MRTI = readVTK_SJF('temperature',161);
%array_uncert = readVTK_SJF('uncertainty',161);

%load 'referenceBased.mat'
%MRTI = (temperature);

% Crop the displayed VOI
MRTI_crop = MRTI;                        % Start with original MRTI
MRTI_crop ( 1 : (VOI.x(1)-1)   , :,:,:) = 0; % Crop top
MRTI_crop ( (VOI.x(2)+1) : end , :,:,:) = 0; % Overwrite MRTI_crop, crop bottom
MRTI_crop ( : ,  1 : (VOI.y(1)-1) ,:,:) = 0; % Crop left
MRTI_crop ( : , (VOI.y(2)+1) : end,:,:) = 0; % Crop right

MRTI_crop_small = MRTI;                      % Start with original MRTI
MRTI_crop_small ( 1 : (VOI_small.x(1)-1)   , :,:,:) = 0; % Crop top
MRTI_crop_small ( (VOI_small.x(2)+1) : end , :,:,:) = 0; % Overwrite MRTI_crop, crop bottom
MRTI_crop_small ( : ,  1 : (VOI_small.y(1)-1) ,:,:) = 0; % Crop left
MRTI_crop_small ( : , (VOI_small.y(2)+1) : end,:,:) = 0; % Crop right

figure(1); imagesc ( MRTI (:,:,time) , [40 80] );
figure(2); imagesc ( MRTI_crop (:,:,time) , [40 80] );
figure(3); imagesc ( MRTI_crop_small (:,:,time) , [40 80] );
center
center_in_meters
VOI
cd ( path22 );
end