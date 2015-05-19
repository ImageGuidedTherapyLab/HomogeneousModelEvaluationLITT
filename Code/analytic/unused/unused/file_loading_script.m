fin1 = fopen('index.txt');

input_cell1 {1,1} = fgetl (fin1); % index
input_cell1 {2,1} = fgetl (fin1); % study
input_cell1 {3,1} = fgetl (fin1); % series
input_cell1 {4,1} = fgetl (fin1); % opttype

fclose(fin1);
path_base = strcat( './workdir/Study',num2str(input_cell1 {2,1}), '/',num2str(input_cell1 {3,1}), '/opt/' );
input_mat = load( strcat( path_base, 'optpp_pds.', input_cell1{4,1}, '.in.1.mat') );

fin2 = fopen( strcat( path_base, 'optpp_pds.', input_cell1{4,1}, '.in.', num2str(input_cell1{1,1}) ) );
ii = 1;
while ~feof(fin2)
    
    input_cell2 {ii,1} = fgetl(fin2);
    ii = ii + 1;
    
end
clear ii



fin3 = fopen( strcat( path_base, 'setup.ini'));
while ~feof(fin3)
    
    fline = fgetl(fin3);
    if strncmp(fline, 'voi =', 5) == 1
        input_cell3 {1,1} = fline;
    end
    if strncmp(fline, 'history =', 9) == 1
        input_cell3 {2,1} = fline;
    end
    
end

fin4 = fopen( strcat( input_mat.mrti, 'temperature.0000.vtk' ));
jj = 1;
for ii = 1:7
    
    fline = fgetl(fin4);
    if ii == 5
        input_cell4 {jj,1} = fline;
        jj = jj +1;
    end
    if ii == 7
        input_cell4 {jj,1} = fline;
        jj = jj +1;
    end
end
clear ii jj




inputdatavars.UID = input_cell1 {2,1};
inputdatavars.powerhistory = input_cell3{2,1};  % Needs further parsing
inputdatavars.opttype = input_cell1 {4,1};
inputdatavars.voi = input_cell3{1,1}; % Needs further parsing
inputdatavars.dimensions = [256;256;1];
%inputdatavars.spacing = 
