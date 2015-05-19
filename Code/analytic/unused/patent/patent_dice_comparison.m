cd /FUS4/data2/sjfahrenholtz/MATLAB/Tests

close all
clear

fit75=load('plot_pattern_bestfit75.mat');
fit76=load('plot_pattern_bestfit76.mat');
fit77=load('plot_pattern_bestfit77.mat');
amalgam1=load('plot_pattern_amalgam1.mat');
amalgam2=load('plot_pattern_amalgam_76_77.mat');

CI_list7 = cell(6,3);
CI_list7{1,1} = 'fit75';
CI_list7{2,1} = 'fit76';
CI_list7{3,1} = 'fit77';
CI_list7{4,1} = 'amalgam1';
CI_list7{5,1} = 'amalgam2_Lo';
CI_list7{6,1} = 'amalgam2_Hi';
CI_list7{1,2} = fit75.hh7.ci(1);
CI_list7{2,2} = fit76.hh7.ci(1);
CI_list7{3,2} = fit77.hh7.ci(1);
CI_list7{4,2} = amalgam1.hh7.ci(1);
CI_list7{5,2} = amalgam2.hhLo7.ci(1);
CI_list7{6,2} = amalgam2.hhHi7.ci(1);

CI_list7{1,3} = fit75.dice_LOOCV7_stats;
CI_list7{2,3} = fit76.dice_LOOCV7_stats;
CI_list7{3,3} = fit77.dice_LOOCV7_stats;
CI_list7{4,3} = amalgam1.dice_LOOCV7_stats;
CI_list7{5,3} = amalgam2.dice_LOOCV7_statsLo;
CI_list7{6,3} = amalgam2.dice_LOOCV7_statsHi;

CI_list3 = cell(6,3);
CI_list3{1,1} = 'fit75';
CI_list3{2,1} = 'fit76';
CI_list3{3,1} = 'fit77';
CI_list3{4,1} = 'amalgam1';
CI_list3{5,1} = 'amalgam2_Lo';
CI_list3{6,1} = 'amalgam2_Hi';
CI_list3{1,2} = fit75.hh3.ci(1);
CI_list3{2,2} = fit76.hh3.ci(1);
CI_list3{3,2} = fit77.hh3.ci(1);
CI_list3{4,2} = amalgam1.hh3.ci(1);
CI_list3{5,2} = amalgam2.hhLo3.ci(1);
CI_list3{6,2} = amalgam2.hhHi3.ci(1);

CI_list3{1,3} = fit75.dice_LOOCV3_stats;
CI_list3{2,3} = fit76.dice_LOOCV3_stats;
CI_list3{3,3} = fit77.dice_LOOCV3_stats;
CI_list3{4,3} = amalgam1.dice_LOOCV3_stats;
CI_list3{5,3} = amalgam2.dice_LOOCV3_statsLo;
CI_list3{6,3} = amalgam2.dice_LOOCV3_statsHi;




