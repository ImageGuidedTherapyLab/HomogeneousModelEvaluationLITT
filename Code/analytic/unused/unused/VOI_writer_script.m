% This is the superscript that allows you to input the VOI information.

% Paths. Be in the folder you want to process. E.g.:
% '/FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/workdir/Patient0002/000/opt'

cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
setenv ( 'WORKDIR' , pwd) ;
work_dir = getenv ( 'WORKDIR' );

patient_index = 1;
Patient_Paths = importdata( 'patient_paths.txt' );
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);

% Patient0002/000
VOI_pre.x = [ 83 123] ;
VOI_pre.y = [ 108 148 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 103 128 0 ];
VOI_pre.time = 60;
VOI = char( strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']') );
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0002/001
VOI_pre.x = [ 89 131 ] ;
VOI_pre.y = [ 105 147 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 110 126 0 ];
VOI_pre.time = 67;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0002/010
VOI_pre.x = [ 100 130 ] ;
VOI_pre.y = [ 110 140 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 115 125 0 ];
VOI_pre.time = 52;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0002/0017
VOI_pre.x = [ 107 137 ] ;
VOI_pre.y = [ 108 139 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 122 123 0 ];
VOI_pre.time = 71;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0002/020
VOI_pre.x = [ 104 134 ] ;
VOI_pre.y = [ 109 140 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 121 124 0 ];
VOI_pre.time = 12;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0002/021
VOI_pre.x = [ 86 108 ] ;
VOI_pre.y = [ 118 140 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 97 129 0 ];
VOI_pre.time = 84;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0002/022
VOI_pre.x = [ 105 135 ] ;
VOI_pre.y = [ 109 139 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 120 124 0 ];
VOI_pre.time = 84;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0003/002
VOI_pre.x = [ 107 137 ] ;
VOI_pre.y = [ 111 141 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 122 126 0 ];
VOI_pre.time = 63;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0003/013
VOI_pre.x = [ 100 140 ] ;
VOI_pre.y = [ 90 120 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 125 106 0 ];
VOI_pre.time = 38;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0003/018
VOI_pre.x = [ 110 140 ] ;
VOI_pre.y = [ 108 132 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 125 118 0 ];
VOI_pre.time = 87;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0003/025
VOI_pre.x = [ 105 139 ] ;
VOI_pre.y = [ 110 142 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 122 126 0 ];
VOI_pre.time = 29;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0004/003
VOI_pre.x = [ 130 170 ] ;
VOI_pre.y = [ 118 148 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 150 133 0 ];
VOI_pre.time = 86;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0004/006
VOI_pre.x = [ 130 160 ] ;
VOI_pre.y = [ 113 135 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 145 120 0 ];
VOI_pre.time = 110;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0005/004
VOI_pre.x = [ 121 145 ] ;
VOI_pre.y = [ 136 160 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 133 143 0 ];
VOI_pre.time = 81;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0005/008
VOI_pre.x = [ 125 145 ] ;
VOI_pre.y = [ 141 163 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 135 156 0 ];
VOI_pre.time = 125;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0005/012
VOI_pre.x = [ 122 142 ] ;
VOI_pre.y = [ 135 153 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 132 143 0 ];
VOI_pre.time = 66;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0006/005 (this patient doesn't seem to have a laser pulse)
VOI_pre.x = [ 100 105 ] ;
VOI_pre.y = [ 100 105 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 102 102 0 ];
VOI_pre.time = 3;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0006/007
VOI_pre.x = [ 100 138 ] ;
VOI_pre.y = [ 130 161 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 119 149 0 ];
VOI_pre.time = 139;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0006/009
VOI_pre.x = [ 112 138 ] ;
VOI_pre.y = [ 150 176 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 126 163 0 ];
VOI_pre.time = 75;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0006/019 (this patient doesn't seem to reach max heating)
VOI_pre.x = [ 110 130 ] ;
VOI_pre.y = [ 137 157 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 120 147 0 ];
VOI_pre.time = 43;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0006/023
VOI_pre.x = [ 104 130 ] ;
VOI_pre.y = [ 150 175 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 117 164 0 ];
VOI_pre.time = 82;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0006/024
VOI_pre.x = [ 90 112 ] ;
VOI_pre.y = [ 137 158 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 99 147 0 ];
VOI_pre.time = 40;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0007/011
VOI_pre.x = [ 88 102 ] ;
VOI_pre.y = [ 70 84 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 96 76 0 ];
VOI_pre.time = 41;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0007/015 (no visible pulse in the MRTI)
VOI_pre.x = [ 100 105 ] ;
VOI_pre.y = [ 100 105 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 102 102 0 ];
VOI_pre.time = 3;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0007/026
VOI_pre.x = [ 95 105 ] ;
VOI_pre.y = [ 70 85 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 100 77 0 ];
VOI_pre.time = 112;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0008/014
VOI_pre.x = [ 91 107 ] ;
VOI_pre.y = [ 131 146 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 99 139 0 ];
VOI_pre.time = 28;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0008/016
VOI_pre.x = [ 105 123 ] ;
VOI_pre.y = [ 129 145 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 133 137 0 ];
VOI_pre.time = 15;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

% Patient0008/027
VOI_pre.x = [ 103 121 ] ;
VOI_pre.y = [ 125 148 ] ;
VOI_pre.z = [0 0];
VOI_pre.center = [ 110 138 0 ];
VOI_pre.time = 110;
patient_index = patient_index + 1;
setenv ( 'PATHPT' , char ( Patient_Paths ( patient_index ) ) );
pathpt = getenv ( 'PATHPT' );
whole_path = strcat ( work_dir, pathpt);
cd (whole_path);
VOI = strcat ( 'voi = [', num2str(VOI_pre.x(1)), ', ', num2str(VOI_pre.x(2)), ', ', num2str(VOI_pre.y(1)), ', ', num2str(VOI_pre.y(2)), ', ', num2str(VOI_pre.z(1)), ', ', num2str(VOI_pre.z(2)), ']');
fid = fopen( 'VOI.txt' ,'wt','b');
fprintf( fid, '%s', VOI);
fclose( fid );

clear