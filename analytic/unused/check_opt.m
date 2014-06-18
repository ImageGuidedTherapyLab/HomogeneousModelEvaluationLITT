function y = check_opt( path_base, opttype )

close all

load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );

% index = load ( 'index.txt' );

[L2norm,dice, tmap_model, MRTI_crop] =  fast_temperature_obj_fxn_sanity ( inputdatavars, 1 );

% index = index + 1;
% csvwrite ('index.txt' , index);
metric(1) = L2norm;
metric(2) = 1 - dice;

%metric = 1 -dice;

% fout = fopen( strcat( file_base,'.out.',num2str(inputdatavars.fileID) ), 'w' );
% fprintf(fout, '%s\n', num2str(metric(1)) );
% fprintf(fout, '%s\n'  , num2str(metric(2)) );
% fclose(fout);

figure(1); imagesc(tmap_model);
figure(2); imagesc(MRTI_crop);

model_deg_threshold = tmap_model >= 57;
MRTI_deg_threshold = MRTI_crop >= 57;
n_model = sum(sum( model_deg_threshold ));
n_MRTI = sum(sum( MRTI_deg_threshold ));
intersection = model_deg_threshold + MRTI_deg_threshold;
intersection2 = intersection > 1;
n_intersection = sum(sum( intersection2 ));
dice = 2*n_intersection / (n_model + n_MRTI) 

figure(3); imagesc( intersection)

y =  metric;