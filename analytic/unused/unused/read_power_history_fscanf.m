% This script is meant to find the powers of the Visualase patients and
% then write them. It should be noted that there is an assumed 0 power at
% 0 time at the beginning of written data.
cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/000/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/000/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/001/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/001/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/010/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/010/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/017/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/017/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/020/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/020/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/021/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/021/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0002/022/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/022/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0003/002/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0003/002/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0003/013/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0003/013/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0003/018/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0003/018/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0003/025/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0003/025/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0004/003/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0004/003/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0004/006/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0004/006/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0005/004/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0005/004/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0005/008/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0005/008/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0005/012/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0005/012/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0006/005/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0006/005/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0006/007/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0006/007/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0006/009/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0006/009/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0006/019/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0006/019/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0006/023/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0006/023/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0006/024/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0006/024/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0007/015/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0007/015/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0007/026/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0007/026/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0008/014/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0008/014/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0008/016/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0008/016/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );

cd '/FUS4/data2/BioTex/BrainNonMDA/processed/Patient0008/027/laser_log'
load 'power_log.txt';
cd '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0008/027/opt'
[Power_intervals] = power_parser_write_DF_array(power_log);
fileID = fopen( 'Power_History.txt' , 'w+' ); fprintf( fileID, Power_intervals); fclose( fileID );
clear