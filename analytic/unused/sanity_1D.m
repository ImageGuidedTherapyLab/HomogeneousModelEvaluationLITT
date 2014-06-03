% This is a sanity script for the 1-D case
clear all
close all
ind_case = 2;
cool_case = 1;
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation

% 1D, mu_eff = 180; cooling u0 = 21 C
opt_type = 'sanity_1D' ;
Study_paths {1,1} = 'Study0035';
Study_paths {1,2} = '0530';
if ind_case ==1;
    Study_paths {1,3} = '1';
elseif ind_case ==2;
    Study_paths {1,3} = '43';
end
param_file = strcat( 'workdir/', Study_paths{1,1}, '/', Study_paths{1,2}, '/opt/', opt_type, '.in.', Study_paths{1,3} );
python_command = strcat( 'unix(''python ./brainsearch.py --param_file ./', param_file, ''')');   % unix(''python test_saveFile.py'')
evalc(python_command);

aa = load( 'TmpDataInput.mat' );
aa.cv.body_temp = str2num( aa.cv.body_temp );
Pwr = 12;
R1 = 0.0007;
R2 = 1;
if cool_case == 1;
    aa.cv.probe_init = str2num( aa.cv.probe_init );
elseif cool_case ==2;
    aa.cv.probe_init = 37;
end

% short range
rr = linspace(0.0001,0.0007,27);
temperature = zeros(length(rr),1);
for ii=1:length(rr);
    [temperature(ii)]= sammodel1D (aa.cv.probe_init , aa.cv.body_temp, aa.cv.k_0, aa.cv.w_0, Pwr, rr(ii), aa.cv.mu_a, aa.cv.mu_s, R1, R2, aa.cv.anfact);
end
figure;plot(rr,temperature);

% Long range
rr = linspace(0.0007,1,9994);
temperature = zeros(length(rr),1);
for ii=1:length(rr);
    [temperature(ii)]= sammodel1D (aa.cv.probe_init , aa.cv.body_temp, aa.cv.k_0, aa.cv.w_0, Pwr, rr(ii), aa.cv.mu_a, aa.cv.mu_s, R1, R2, aa.cv.anfact);
end
figure;plot(rr,temperature); % Whole space
figure;plot(rr(1:1000),temperature(1:1000)); % first 10 cm = 0.1 m
figure;plot(rr(1:400),temperature(1:400));
index = find(temperature(100:end)<57,1,'first') + 99;
rr(index)
max(temperature)