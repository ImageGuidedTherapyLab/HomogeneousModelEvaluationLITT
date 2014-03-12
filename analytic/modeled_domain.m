%This function accepts the FOV, matrix, scaling, and the number of image pixels that will be modeled and returns a
%suggested domain for the Bioheat1D.m function. The FOV is a structure
%describing the number of x and y pixels and the number of slices. The the output domain has
%tighter spacing as the scaling increases. I.e. the number of points is
%linear with the scaling. Scaling is a 3-field structure too (for x,y,z).
%'mod_point' is the number of image pixels that are meant to be modeled.

function [dom,MRTI_pix,mod_pix]=modeled_domain(FOV,matrix,scaling,mod_point)

%Modeled domain in meters
dom.x=FOV.x*mod_point.x/matrix.x;  %dom = FOV * (fraction modeled) = FOV *('modeled points'/'image matrix points')
dom.y=FOV.y*mod_point.y/matrix.y;
dom.z=FOV.z*mod_point.z/matrix.z;

%Modeled domain points
dom.pointy=floor(mod_point.x*scaling.x);
dom.pointx=floor(mod_point.y*scaling.y);
dom.pointz=floor(mod_point.z*scaling.z);

%Find the MRTI pixel size
MRTI_pix.x=FOV.x/matrix.x;
MRTI_pix.y=FOV.y/matrix.y;
MRTI_pix.z=FOV.z/matrix.z;

%Modeled pixel size
mod_pix.x=dom.x/dom.pointx;
mod_pix.y=dom.y/dom.pointy;
mod_pix.z=dom.z/dom.pointz;

end