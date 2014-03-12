cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

setenv ( 'PATHPT' , '/workdir/Patient0006/007/opt' );
pathpt = getenv ( 'PATHPT' );

load index.txt

[metric] = obj_fxn_unified_optics ( path22, pathpt, index );

index = index + 1;
csvwrite ('index.txt' , index);