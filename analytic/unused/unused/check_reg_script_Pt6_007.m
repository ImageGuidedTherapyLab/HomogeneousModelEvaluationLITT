% This is the superscript that has the paths for all of the patient data.

% Paths. Be in the folder you want to process. E.g.:
% '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/000/opt'

cd /FUS4/data2/BioTex/BrainNonMDA/processed/Patient0006/007/laser
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

time = 143;
crop_big = 25;
crop_small = 5;

[ center, VOI ] = check_reg ( path22 ,time,crop_big,crop_small );