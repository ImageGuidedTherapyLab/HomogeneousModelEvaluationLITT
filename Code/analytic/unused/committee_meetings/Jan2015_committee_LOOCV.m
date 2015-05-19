close all
clear
clc

cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
load Jan2015_committee_all_27.mat

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

% LOOCV, bad;
% [val index] = min(LOOCV.dice.values{1});
% [val index] = max(opt.mu_eff.values{1});
index = 2;
[nm1 nm2] = LOOCV.paths.paths{1}{index,:};
Study_paths = {nm1 nm2};
mu_eff_list = opt.mu_eff.values{1};
mu_eff_list(index) = [];
mu_eff = num2str( mean(mu_eff_list) );

% % Opt, bad;
% index = 2;
% [nm1 nm2] = LOOCV.paths.paths{1}{index,:};
% Study_paths = {nm1 nm2};
% mu_eff = num2str( opt.mu_eff.values{1}(index) );

[L2norm, dice, tmap_model, MRTI_crop,mat_struct,quality] = check_opt_GoodBadUgly22( Study_paths, opttype, mu_eff);

figure(1); imagesc(tmap_model{1}, [35 85]);
figure(2); imagesc(MRTI_crop{1}, [35 85]);

sz = size(tmap_model{1})
sz(1).*mat_struct.file.inputdatavars.spacing(1)
sz(2).*mat_struct.file.inputdatavars.spacing(2)

% Bad Optimization pic

dice
Study_paths
opt.mu_eff.values{1}(index)
index

model_deg_threshold = (tmap_model{1} >= 57);
MRTI_deg_threshold = (MRTI_crop{1} >= 57) .* 2;
merge = model_deg_threshold + MRTI_deg_threshold;
figure(3); imagesc(merge, [0 3]);