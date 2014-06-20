function y = check_opt( paths, opttype, best_iter )

num_studies = size(paths);



for ii = 1:num_studies
    clear inputdatavars
    path_base = strcat ('./workdir/', paths{ii,1},'/',paths{ii,2},'/opt/optpp_pds.',opttype,'.in.', num2str(best_iter),'.mat');
    
    load ( path_base );

    



%load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );

% index = load ( 'index.txt' );

inputdatavars.cv.mu_eff_healthy = num2str( mu_eff );

[L2norm(ii,:) , dice(ii,:), tmap_model(:,:,ii), MRTI_crop(:,:,ii)] =  fast_temperature_obj_fxn_sanity ( inputdatavars, 1 );

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