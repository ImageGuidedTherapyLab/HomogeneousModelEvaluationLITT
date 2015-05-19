% plots for paper

temp_paths=Study_paths25;
toss_indexLOOCV=find(dice_LOOCV25<.7);
temp_paths(toss_indexLOOCV,:)=[];

best_iter_pass_LOOCV25 = best_iter25;
best_iter_pass_LOOCV25(toss_indexLOOCV,:)=[];
[L2norm0455, dice0455, tmap_model0455, MRTI_crop0455,mat_struct0455,quality0455] = check_opt_GoodBadUgly( {temp_paths{5,1}, temp_paths{5,2}}, opttype, best_iter_pass_LOOCV25(5) );

model_deg_threshold = tmap_model0455{1} >= 57;
MRTI_deg_threshold = 2*(MRTI_crop0455{1} >= 57);
n_model = sum(sum( model_deg_threshold ));
n_MRTI = sum(sum( MRTI_deg_threshold ));
sum_of_overlap = model_deg_threshold + MRTI_deg_threshold;
intersection = sum_of_overlap > 1;
n_intersection = sum(sum( intersection ));
dice = 2*n_intersection / (n_model + n_MRTI) ;