function [L2norm, dice, tmap_model, MRTI_crop,mat_struct] = check_opt_LOOCV( paths, opttype, mu_eff_LOOCV )

num_studies = size(paths,1);
L2norm = zeros(num_studies,1);
dice = zeros(num_studies,1);
tmap_model = cell(num_studies,1);
MRTI_crop = cell(num_studies,1);
mat_struct( num_studies ) = struct ();

for ii = 1:num_studies
    clear inputdatavars
    path_base = strcat ('./workdir/', paths{ii,1},'/',paths{ii,2},'/opt/optpp_pds.',opttype,'.in.1.mat');
    
    load ( path_base );
    mu_eff_iter = mu_eff_LOOCV;
    mu_eff_iter ( ii ) = []; % Remove the 0
    mu_eff_iter = mean ( mu_eff_iter );
    inputdatavars.cv.mu_eff_healthy = num2str( mu_eff_iter);

    [L2norm(ii) , dice(ii), tmap_model{ii}, MRTI_crop{ii}] =  fast_temperature_obj_fxn_sanity ( inputdatavars, 1 );
    
end
end

% figure(1); imagesc(tmap_model);
% figure(2); imagesc(MRTI_crop);
% 
% model_deg_threshold = tmap_model >= 57;
% MRTI_deg_threshold = MRTI_crop >= 57;
% n_model = sum(sum( model_deg_threshold ));
% n_MRTI = sum(sum( MRTI_deg_threshold ));
% intersection = model_deg_threshold + MRTI_deg_threshold;
% intersection2 = intersection > 1;
% n_intersection = sum(sum( intersection2 ));
% dice = 2*n_intersection / (n_model + n_MRTI) 
% 
% figure(3); imagesc( intersection)

