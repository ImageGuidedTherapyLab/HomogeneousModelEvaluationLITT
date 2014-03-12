% This is the superscript that has the paths for all of the patient data.
% Set the index to whatever the best iteration was.
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

setenv ( 'PATHPT' , '/workdir/Patient0006/007/opt' );
pathpt = getenv ( 'PATHPT' );

path_total = strcat ( path22, pathpt );
cd (path_total)
load ('obj_fxn_list.txt');
index = min ( obj_fxn_list );
cd ( path22 );

[metric,matched_mod,diff_Iso] = unified_optics_recover_VTK ( path22, pathpt, index );