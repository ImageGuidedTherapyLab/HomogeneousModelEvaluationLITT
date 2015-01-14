%This function accepts the FOV, matrix, scaling, and the number of image pixels that will be modeled and returns a
%suggested domain for the Bioheat1D.m function. The FOV is an array
%describing the length of x and y dimensions and the number of slices. The output domain has
%tighter spacing as the scaling increases. I.e. the number of points is
%linear with the scaling. Scaling is a 3-field array too (for x,y,z).
%'mod_point' is the number of image pixels that are meant to be modeled.

function [dom,dom_point,MRTI_pix,mod_pix]=modeled_domain_array(FOV,matrix,scaling,mod_point)

%Modeled domain in meters
dom=FOV.*mod_point./matrix;  %dom = FOV * (fraction modeled) = FOV *('modeled points'/'image matrix points')

%Modeled domain points
dom_point=mod_point.*scaling;

%Find the MRTI pixel size
MRTI_pix=FOV./matrix;

%Modeled pixel size
mod_pix=dom./dom_point;

end