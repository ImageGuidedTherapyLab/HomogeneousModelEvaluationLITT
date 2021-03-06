% y = f(x) for Rosenbrock
function y = f(~)
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation/
% setenv ( 'PATH22' , pwd);
% path22 = getenv ( 'PATH22' );

%file_loading_script

% save('./superTemp.mat','x');
% disp(x)

inputdatavars = load('./TmpDataInput.mat');

% index = load ( 'index.txt' );

%[L2norm,dice, ~,~] =  fast_temperature_obj_fxn_sanity ( inputdatavars, 1 );
[L2norm,dice, ~,~] =  temperature_obj_fxn ( inputdatavars, 10 );
% index = index + 1;
% csvwrite ('index.txt' , index);
metric(1) = L2norm;
%metric(2) = 1 - dice;
%metric(isnan(metric)) = 0;
metric(2) = 1/(dice + 1E-7);
metric(3) = dice;
metric(4) = -dice;
metric(5) = nthroot(metric(1),4)/500 + metric(4).*500;
%metric(5) = metric(1) + metric(2);
%metric = 1 -dice;

file_base = strcat( './workdir/',inputdatavars.patientID,'/',inputdatavars.UID,'/opt/optpp_pds.',inputdatavars.opttype);
fout = fopen( strcat( file_base,'.out.',num2str(inputdatavars.fileID) ), 'w' );
fprintf(fout, '%s\n', num2str(metric(1)) );
fprintf(fout, '%s\n', num2str(metric(2)) );
fprintf(fout, '%s\n', num2str(metric(3)) );
fprintf(fout, '%s\n', num2str(metric(4)) );
fprintf(fout, '%s\n', num2str(metric(5)) );
fclose(fout);

save ( strcat( file_base,'.in.',num2str(inputdatavars.fileID), '.mat'), 'inputdatavars');



y =  metric;


end
% patientID = strcat ( inputdatavars.patientID, '/', inputdatavars.patientID, '/');
% patient_index = load ( 'patient_index.txt' );
% vtk_times = load ( 'VTK_patient_times.txt' );
% 
% Patient_Paths = importdata( 'patient_paths.txt' );
% 
% setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
% pathpt = getenv ( 'PATHPT' );
% 
% vtk_index = vtk_times ( patient_index );

% [metric] =  fast_temperature_obj_fxn ( path22, pathpt, index, vtk_index );
% x(1)
% disp('k_0')
% x(2)
% disp('mu_a')
% x(3)