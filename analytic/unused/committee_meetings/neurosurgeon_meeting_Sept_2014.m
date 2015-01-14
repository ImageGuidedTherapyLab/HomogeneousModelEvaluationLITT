close all
clear

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

aaa = load ('./workdir/Study0026/0450/opt/optpp_pds.bestfit76.in.10.mat');

bbb = cell(20,1);
L2norm = zeros(20,1);
dice = zeros(20,1);
tmap_model = cell(20,1);
MRTI_crop = cell(20,1);
mu_eff_samples = linspace(100,2000,20);



for ii = 1:length(bbb)
    bbb{ii} = aaa;
    bbb{ii}.inputdatavars.cv.mu_eff_healthy = num2str(mu_eff_samples(ii));
    
    [L2norm(ii), dice(ii), tmap_model{ii}, MRTI_crop{ii}] = temperature_obj_fxn ( bbb{ii}.inputdatavars, 10 );
    
    figure; imagesc(tmap_model{ii},[30 80]);
    
end


% Run 2nd
ccc = load ('./workdir/Study0026/0450/opt/optpp_pds.bestfit76.in.15.mat');

[L2norm, dice, tmap_model, MRTI_crop] = temperature_obj_fxn ( ccc.inputdatavars, 10 );

figure; imagesc(tmap_model,[30 80]); colorbar;
figure; imagesc(MRTI_crop, [30 80]); colorbar;

model_deg_threshold = tmap_model >= 57;
MRTI_deg_threshold = MRTI_crop >= 57;
n_model = sum(sum( model_deg_threshold ));
n_MRTI = sum(sum( MRTI_deg_threshold ));
intersection = model_deg_threshold + 2*MRTI_deg_threshold;

figure; imagesc(intersection,[0,3]); colorbar



    