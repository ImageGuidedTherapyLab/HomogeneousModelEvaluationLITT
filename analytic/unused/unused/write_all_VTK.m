cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

Patient_Paths = importdata( 'patient_paths.txt' ); % Load the paths to the datasets
NumPatients = size ( Patient_Paths, 1 );  % Read in the total number of datasets
patient_index = 1; % Initialize the patient path index

for ii = 1 : NumPatients  % Iterate over all patient datasets
        
    setenv ( 'PATHPT' , char ( Patient_Paths ( ii ) ) ); % Set the patient path using the patient path index
    pathpt = getenv ( 'PATHPT' );
    
    path_total = strcat ( path22, pathpt );
    cd (path_total)
    load ('obj_fxn_list.txt'); % Read in all of the output objective functions
    [~,index] = min ( obj_fxn_list ); % Identify the best iteration.
    
    [~] = fast_temperature_Write_all ( path22, pathpt, index );
    
end