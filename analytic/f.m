% y = f(x) for Rosenbrock
function y = f(~)
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation/
% setenv ( 'PATH22' , pwd);
% path22 = getenv ( 'PATH22' );

inputdatavars = load('./TmpDataInput.mat');

% index = load ( 'index.txt' );

[metric,~,~] =  fast_temperature_obj_fxn ( inputdatavars );

% index = index + 1;
% csvwrite ('index.txt' , index);

y =  metric;


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