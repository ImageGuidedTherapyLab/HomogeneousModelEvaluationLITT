close all
clear
clc

cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
load Jan2015_committee.mat

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

% Good optimization pic
Study_paths={'Study0030' '0495'};
mu_eff = num2str(opt.mu_eff.all.values(1));

% % Bad Optimization pic
% Study_paths={'Study0021' '0415'};
% mu_eff = num2str(opt.mu_eff.all.values(22));

[L2norm, dice, tmap_model, MRTI_crop,mat_struct,quality] = check_opt_GoodBadUgly22( Study_paths, opttype, mu_eff);

figure(1); imagesc(tmap_model{1}, [35 85]);
figure(2); imagesc(MRTI_crop{1}, [35 85]);

sz = size(tmap_model{1})
sz(1).*mat_struct.file.inputdatavars.spacing(1)
sz(2).*mat_struct.file.inputdatavars.spacing(2)

% Bad Optimization pic

dice

model_deg_threshold = (tmap_model{1} >= 57);
MRTI_deg_threshold = (MRTI_crop{1} >= 57) .* 2;
merge = model_deg_threshold + MRTI_deg_threshold;
figure(3); imagesc(merge);


% Histogram, wider;
bins = (1:10).*40;
bins = bins - 20;
figure(4);hist(opt.mu_eff.values{1},bins);

% Histogram, small;
figure(5);hist(opt.mu_eff.values{1});


% LOOCV, bad;
[val index] = min(LOOCV.dice.values{1});
Study_paths = {'Study0026' '0447'};
mu_eff_list = opt.mu_eff.values{1};
mu_eff_list(index) = [];
mu_eff = num2str( mean(mu_eff_list) );
[L2norm, dice, tmap_model, MRTI_crop,mat_struct,quality] = check_opt_GoodBadUgly22( Study_paths, opttype, mu_eff);


%[tmap, MRTI_crop] = opt_check_Fxn( Study_paths, opttype, opt);