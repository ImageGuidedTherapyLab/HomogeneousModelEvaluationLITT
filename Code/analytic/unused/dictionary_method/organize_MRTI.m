% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [MRTI_crop, MRTI_isotherm, MRTI_list, n_MRTI] = organize_MRTI ( inputdatavars );
% Record the working directory
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% Make the path to the patient directory
patientID = strcat ( inputdatavars.patientID, '/', inputdatavars.UID, '/');
patient_opt_path = strcat ( path22, '/workdir/', patientID, 'opt/' ); % The '/workdir/' needs the first backslash coz it uses absolute path

% Make the path to the VTK
patient_MRTI_path = strcat ( 'StudyDatabase/', patientID, 'vtk/referenceBased/' );

% Change directory and load the temperature from VTK
cd (patient_MRTI_path);


% Make the VOI; Note the 'inputdatavars.voi' is from ParaView.
VOI.x = double( inputdatavars.voi(3:4)); % The weird index assignment is coz it's from ParaView.
VOI.y = double( inputdatavars.voi(1:2));
VOI.z = double( inputdatavars.voi(5:6));

if sum( inputdatavars.UID == '0496' ) ==4 %
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0402' ) ==4%
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 3;
elseif sum( inputdatavars.UID == '0389' ) ==4%
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 2;
elseif sum( inputdatavars.UID == '0385' ) ==4%
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y - 1;
    inputdatavars.maxheatid = 109;
elseif sum( inputdatavars.UID == '0476' ) ==4%
    VOI.x = VOI.x + 5;
    VOI.y = VOI.y - 5;
elseif sum( inputdatavars.UID == '0477' ) ==4%
    VOI.x = VOI.x + 3;
    VOI.y = VOI.y - 2;
elseif sum( inputdatavars.UID == '0438' ) ==4%
    VOI.x = VOI.x - 1;
    VOI.y = VOI.y + 0;
elseif sum( inputdatavars.UID == '0435' ) ==4%
    VOI.x = VOI.x - 3;
    VOI.y = VOI.y + 3;
elseif sum( inputdatavars.UID == '0436' ) ==4%
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 6;
    inputdatavars.maxheatid = 39;
elseif sum( inputdatavars.UID == '0466' ) ==4%
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0468' ) ==4%
    VOI.x = VOI.x - 3;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0471' ) ==4%
    VOI.x = VOI.x + 3;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0447' ) ==4%
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0453' ) ==4%
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0451' ) ==4%
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0418' ) ==4%
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 0;
    %inputdatavars.maxheatid = 85;
elseif sum( inputdatavars.UID == '0409' ) ==4%
    VOI.x = VOI.x - 2;
    VOI.y = VOI.y + 11;
elseif sum( inputdatavars.UID == '0414' ) ==4%
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0415' ) ==4%
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 1;
end

MRTI = readVTK_one_time('temperature', inputdatavars.maxheatid);   % This 'vtkNumber' should b
cd (patient_opt_path);

% Crop the MRTI using the VOI.
% This is the VOI.x and VOI.y swapped for ParaView
MRTI ( 1:(VOI.y(1)-1), :, : ) = 0;  % The -1 and +1 make sure the VOI indices aren't cut
MRTI ( (VOI.y(2)+1):end, :, : )  = 0;
MRTI ( :, 1:(VOI.x(1)-1), : ) = 0;
MRTI ( :,(VOI.x(2)+1):end, : ) = 0;

MRTI_crop = MRTI( (VOI.y(1) ):(VOI.y(2) ) , (VOI.x(1) ):(VOI.x(2) ) ); % Set the cropped region
%MRTI_crop = permute( MRTI_crop, [2 1 3]);

cd (path22);

MRTI_isotherm = zeros( size(MRTI_crop,1), size(MRTI_crop,2) ,15);
n_MRTI = zeros(15,1);
MRTI_list = cell(15,1);

for kk=1:15
    MRTI_tmp = MRTI_crop >= ( 50 + kk);
    MRTI_isotherm(:,:,kk) = imfill( ExtractNLargestBlobs(MRTI_tmp,1) , 'holes')  ;  % 1st, it keeps the largest contiguous region, then it fills in any holes
    [MRTI_row, MRTI_column] = find( MRTI_isotherm(:,:,kk) ==1);
    MRTI_list{kk} = [(MRTI_row .* inputdatavars.spacing(1)) (MRTI_column .* inputdatavars.spacing(2))];
    n_MRTI(kk) =  sum( sum( MRTI_isotherm(:,:,kk)  ));
end
clear kk

end