function [dice, map1_threshold, map2_threshold, label] = dice_maps (input1, input2, isotherm);

map1_threshold = input1 >= isotherm;
map2_threshold = input2 >= isotherm;
n_map1 = sum(sum( map1_threshold ));
n_map2 = sum(sum( map2_threshold ));
intersection = map1_threshold + map2_threshold;
intersection = intersection > 1;
n_intersection = sum(sum( intersection ));
dice = 2*n_intersection / (n_map1 + n_map2) ;

label = map1_threshold + 2.*map2_threshold;
end