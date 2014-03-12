% This is the superscript that has the paths for all of the patient data.

% Paths. Be in the folder you want to process. E.g.:
% '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/000/opt'
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

setenv ( 'PATHPT' , '/workdir/Patient0006/007/opt' );
pathpt = getenv ( 'PATHPT' );

%load index.txt
index = 1;
[metric] = fast_temperature_obj_fxn ( path22, pathpt, index );

% index = index + 1;
% csvwrite ('index.txt' , index);