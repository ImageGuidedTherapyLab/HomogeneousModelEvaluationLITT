% script to compare FEM and GF kernels, using Study0030/0488
close all
clear

cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests
load 'SJF_manuscript_fishing.mat';
%DF_data = load('forsam.mat');
DF_data = load('DF_PlanningWorkspace14Nov2014.mat');
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

% survival_plot22 (opt.dice.all.values, opt.dice.all, fig_labels, LOOCV.run1, LOOCV.run2);
% 
% survival_plot22 (DF_data.opt.dice.all.values, DF_data.opt.dice.all, DF_data.fig_labels, DF_data.LOOCV.run1, DF_data.LOOCV.run2);

opt_check22

GF_threshold = tmap_model{1} >= 57;
MRTI_deg_threshold = MRTI_crop{1} >= 57;

cd /tmp/outputs/dakota/0488
FEM_model = readVTK_one_time22('roisemdose.heating', mat_struct.file.inputdatavars.maxheatid);
FEM_model (FEM_model < 1) = 0;
FEM_model (FEM_model >=1) = 1;

MRTI_dose = readVTK_one_time22('roimrtidose.heating', mat_struct.file.inputdatavars.maxheatid);
MRTI_dose (MRTI_dose < 1) = 0;
MRTI_dose (MRTI_dose >=1) = 1;


survival_plot22 (DF_data.opt.dice.values{1}, opt.dice.values{1}, DF_data.LOOCV.dice.values{1}, LOOCV.dice.values{1});

figure; imagesc(GF_threshold);
figure; imagesc(FEM_model);
%figure; imagesc(MRTI_dose);
figure; imagesc(MRTI_deg_threshold);

% Crappy inverse percentiles

% FEM LOOCV
invprctile(DF_data.LOOCV.dice.values{1},[0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.85])

% FEM opt
invprctile(DF_data.opt.dice.values{1},[0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.85])

%GF LOOCV
invprctile(LOOCV.dice.values{1},[0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.85])

%GF opt
invprctile(opt.dice.values{1},[0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.85])


% intersection = GF_threshold + 2.*FEM_model + 4.*MRTI_deg_threshold;
% figure; imagesc(intersection);
% disp('Dark blue = 0 = all no overlap');
% disp('Mid blue = 1 = Only GF kernel');
% disp('Dark Cyan = 2 = Only FEM kernel');
% disp('Cyan = 3 = GF and FEM kernels');
% disp('Yellow = 4 = MRTI only');
% disp('Bright red = 6 = MRTI + FEM');
% disp('Dark red = 7 = all');

intersection = GF_threshold + 2.*MRTI_deg_threshold;;
figure; imagesc(intersection);
disp('Dark blue = 0 = all no overlap');
disp('Cyan = 1 = Only GF kernel');
disp('Gold = 2 = Only MRTI');
disp('Dark Red = 3 = Both');

intersection = FEM_model + 2.*MRTI_deg_threshold;;
figure; imagesc(intersection);
disp('Dark blue = 0 = all no overlap');
disp('Cyan = 1 = Only FEM kernel');
disp('Gold = 2 = Only MRTI');
disp('Dark Red = 3 = Both');

%0.038417 x 0.033732
