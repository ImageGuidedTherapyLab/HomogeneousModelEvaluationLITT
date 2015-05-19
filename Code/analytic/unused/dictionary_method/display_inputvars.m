function [everything, rotate_size,max_phys_sz] =  display_inputvars ( patient_ix);

% Identify the studies to be examined.
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file

opttype = 'bestfit50' ;

datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);

for ii = 1:num_studies
    
    Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end

clear ii

spacing_holder = zeros( length(patient_ix), 7);
%rot_holder = spacing_holder;
rot_holder = zeros( length(patient_ix), 4);
everything = zeros( length(patient_ix), 9);

for ii = patient_ix
    path_base = strcat ( 'workdir/',Study_paths{ii,1}, '/', Study_paths{ii,2}, '/opt');
    load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );
    spacing_holder(ii,1) = str2num( inputdatavars.UID );
    spacing_holder(ii,2:4)= inputdatavars.spacing;
    spacing_holder(ii,6)  = inputdatavars.voi(2)-inputdatavars.voi(1) + 1;
    spacing_holder(ii,7)  = inputdatavars.voi(4)-inputdatavars.voi(3) + 1;
    
    rot_holder(ii,1) = spacing_holder(ii,1);
    rot_holder(ii,2) = str2num(inputdatavars.cv.x_rotate);
    rot_holder(ii,3) = str2num(inputdatavars.cv.y_rotate);
    rot_holder(ii,4) = str2num(inputdatavars.cv.z_rotate);
    
    everything(ii,1) = spacing_holder(ii,1);
    everything(ii,4:6) = spacing_holder(ii,2:4);
    everything(ii,7:8) = spacing_holder(ii,6:7);
    
end
clear ii

x_unique = unique(spacing_holder(:,2));
x_label  = zeros( length(patient_ix),1);
y_unique = unique(spacing_holder(:,3));
y_label = x_label;
z_unique = unique(spacing_holder(:,4));
z_label = x_label;

for ii = 1:length(x_unique)
    ix =  spacing_holder(:,2) == x_unique(ii);
    x_label( ix) = ii;
end
clear ii
for ii = 1:length(y_unique)
    ix =  spacing_holder(:,3) == y_unique(ii);
    y_label( ix) = ii;
end
clear ii
for ii = 1:length(z_unique)
    ix =  spacing_holder(:,4) == z_unique(ii);
    z_label( ix) = ii;
end
clear ii

combo_label = str2num( strcat( num2str( x_label ), num2str( y_label ), num2str( z_label) ));
combo_unique = unique( combo_label);

for ii = 1:length(combo_unique)
    ix = combo_label == combo_unique(ii);
    spacing_holder(ix,5) = ii;
end

combo_unique = unique( spacing_holder(:,5) );
matrix_size = zeros( length( combo_unique), 3);
matrix_size(:,1) = combo_unique;
for ii = 1:length(combo_unique)
    ix = find( spacing_holder(:,5) == ii);
    matrix_size(ii,2) = max ( spacing_holder( ix , 6));
    matrix_size(ii,3) = max ( spacing_holder( ix , 7));
end


% Now rotations
neg_rot = rot_holder(:,4) < 0;
rot_holder( neg_rot,4) = 360 - abs(rot_holder( neg_rot,4) );
theta_high = rot_holder(:,4) >= 180;
rot_holder( theta_high,4) = rot_holder( theta_high,4) - 180;
everything(:,9) = rot_holder(:,4);

x_unique = unique(rot_holder(:,2));
x_label  = zeros( length(patient_ix),1);
y_unique = unique(rot_holder(:,3));
y_label = x_label;
z_unique = unique(rot_holder(:,4));
z_label = x_label;

% for ii = 1:length(x_unique)
%     ix =  rot_holder(:,2) == x_unique(ii);
%     x_label( ix) = ii;
% end
% clear ii
% for ii = 1:length(y_unique)
%     ix =  rot_holder(:,3) == y_unique(ii);
%     y_label( ix) = ii;
% end
%clear ii
for ii = 1:length(z_unique)
    ix =  rot_holder(:,4) == z_unique(ii);
    z_label( ix) = ii;
end
clear ii

%combo_label = str2num( strcat( num2str( x_label ), num2str( y_label ), num2str( z_label) ));
combo_unique = unique( z_label);

for ii = 1:length(combo_unique)
    ix = z_label == combo_unique(ii);
    rot_holder(ix,5) = ii;
end
clear ii

everything(:,2) =  spacing_holder(:,5) ;
everything(:,3) =  rot_holder(:,5);

num_matrix = max( everything(:,2) );
rotate_size = zeros( num_matrix,3);
rotate_size(:,1) = linspace(1,num_matrix,num_matrix);
for ii = 1:num_matrix
    ix = find( everything(:,2) == ii);
    tmp_sz = zeros( length(ix) , 2);
    for kk = 1:length(ix)
        aa = ones( everything( ix(kk), 7), everything( ix(kk),8));
        bb = imrotate( aa, everything( ix(kk),9));
        tmp_sz(kk,:) = size (bb);
    end
    rotate_size(ii,2) = max(tmp_sz(:,1));
    rotate_size(ii,3) = max(tmp_sz(:,2));
end
clear ii kk

% lx = 0.000001;
% ly = lx;
physical_sz = zeros( length(patient_ix),2);
for ii = patient_ix
    aa = ones( everything(ii, 7), everything(ii,8));
    bb = imrotate(aa, everything( ii,9));
    bb_sz = size(bb);
    physical_sz(ii,1) = bb_sz(1).* everything(ii,4);
    physical_sz(ii,2) = bb_sz(2).* everything(ii,5);

end
max_phys_sz =zeros(2,4);
max_phys_sz(1,1) = max(physical_sz(:,1));
max_phys_sz(2,1) = max(physical_sz(:,2));
max_phys_sz(1,2) = max(everything(:,4));
max_phys_sz(2,2) = max(everything(:,5));
max_phys_sz(1,3) = min(everything(:,4));
max_phys_sz(2,3) = min(everything(:,5));
max_phys_sz(1,4) = max( rotate_size(:,2) );
max_phys_sz(2,4) = max( rotate_size(:,3) );

end