%This function accepts the model prediction and the MRTI and uses user
%measured MRTI heat centroid and the angle that the laser fiber uses.
%Center, start, and stop should be listed as [y x], because of the MATLAB
%column-major convention.

function [rot_model] = rotate_model_data [model_D, MRTI_D, center, start, stop,scaling];

diff = stop - start;
angle=atan(diff(1)/diff(2));
rotation_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)];  %Rotation matrix to make model match MRTI orientation

%Make a large matrix that won't crop the rotated model data
dummy_matrix=zeros((scaling.x*size(MRTI_D,2)),scaling.y*size(MRTI_D,1)); %Reverse the order in order to account for column-major
dummy_center=scaling*center; %Makes the center of the dummy matrix match the MRTI data
model_center=[floor(size(model_D,2)/2) floor(size(model_D/1)/2)]; %Find the center of the model data image

%Write the model data to the dummy matrix
dummy_matrix=