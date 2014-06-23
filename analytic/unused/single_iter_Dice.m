function [ dice, model_deg_threshold, MRTI_deg_threshold, intersection ] = single_iter_Dice (tmap_model, MRTI_crop);

model_deg_threshold = tmap_model >= 57;
MRTI_deg_threshold = MRTI_crop >= 57;
n_model = sum(sum( model_deg_threshold ));
n_MRTI = sum(sum( MRTI_deg_threshold ));

MRTI_label = MRTI_deg_threshold .* 2;

intersection = model_deg_threshold + MRTI_deg_threshold;
intersection = intersection > 1;
intersection_label = model_deg_threshold + MRTI_label;
n_intersection = sum(sum( intersection ));
dice = 2*n_intersection / (n_model + n_MRTI);
figure(4); imagesc(model_deg_threshold);
figure(5); imagesc(MRTI_deg_threshold);
figure(6); imagesc(intersection_label);
figure(7); imagesc(tmap_model,[30 90]);
figure(8); imagesc(MRTI_crop, [30 90]);