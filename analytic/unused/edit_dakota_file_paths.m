clear
close all

Study_paths {1,1} = 'Study0035';
Study_paths {1,2} = '0530';
Study_paths {2,1} = 'Study0023';
Study_paths {2,2} = '0433';
Study_paths {3,1} = 'Study0023';
Study_paths {3,2} = '0428';
Study_paths {4,1} = 'Study0030';
Study_paths {4,2} = '0495';
Study_paths {5,1} = 'Study0030';
Study_paths {5,2} = '0497';
Study_paths {6,1} = 'Study0030';
Study_paths {6,2} = '0488';
Study_paths {7,1} = 'Study0030';
Study_paths {7,2} = '0491';
Study_paths {8,1} = 'Study0030';
Study_paths {8,2} = '0496';
Study_paths {9,1} = 'Study0030';
Study_paths {9,2} = '0490';
Study_paths {10,1} = 'Study0017';
Study_paths {10,2} = '0377';
Study_paths {11,1} = 'Study0017';
Study_paths {11,2} = '0378';
Study_paths {12,1} = 'Study0018';
Study_paths {12,2} = '0402';
Study_paths {13,1} = 'Study0018';
Study_paths {13,2} = '0389';
Study_paths {14,1} = 'Study0018';
Study_paths {14,2} = '0385';
Study_paths {15,1} = 'Study0029';
Study_paths {15,2} = '0476';
Study_paths {16,1} = 'Study0029';
Study_paths {16,2} = '0477';
Study_paths {17,1} = 'Study0025';
Study_paths {17,2} = '0438';
Study_paths {18,1} = 'Study0025';
Study_paths {18,2} = '0435';
Study_paths {19,1} = 'Study0025';
Study_paths {19,2} = '0440';
Study_paths {20,1} = 'Study0025';
Study_paths {20,2} = '0436';
Study_paths {21,1} = 'Study0028';
Study_paths {21,2} = '0466';
Study_paths {22,1} = 'Study0028';
Study_paths {22,2} = '0468';
Study_paths {23,1} = 'Study0028';
Study_paths {23,2} = '0471';
Study_paths {24,1} = 'Study0052';
Study_paths {24,2} = '0725';
Study_paths {25,1} = 'Study0052';
Study_paths {25,2} = '0720';
Study_paths {26,1} = 'Study0026';
Study_paths {26,2} = '0447';
Study_paths {27,1} = 'Study0026';
Study_paths {27,2} = '0457';
Study_paths {28,1} = 'Study0026';
Study_paths {28,2} = '0455';
Study_paths {29,1} = 'Study0026';
Study_paths {29,2} = '0453';
Study_paths {30,1} = 'Study0026';
Study_paths {30,2} = '0450';
Study_paths {31,1} = 'Study0026';
Study_paths {31,2} = '0451';
Study_paths {32,1} = 'Study0022';
Study_paths {32,2} = '0418';
Study_paths {33,1} = 'Study0022';
Study_paths {33,2} = '0417';
Study_paths {34,1} = 'Study0021';
Study_paths {34,2} = '0409';
Study_paths {35,1} = 'Study0021';
Study_paths {35,2} = '0414';
Study_paths {36,1} = 'Study0021';
Study_paths {36,2} = '0415';
Study_paths {37,1} = 'Study0006';
Study_paths {37,2} = '0183';
num_studies = size(Study_paths,1);

dakota_filename_in = 'dakota_bestfit1';
dakota_filename_out= 'dakota_bestfit1';
%file_base = strcat( './workdir/',Study_paths{1,1},'/',Study_paths{1,2},'/opt/');
%copyfile( strcat( file_base,dakota_filename_in,'.in' ), strcat( file_base, dakota_filename_out, '.txt') );   % unix(''python test_sa
% shell_command = strcat( 'unix(''cp workdir/',Study_paths{1,1},'/',Study_paths{1,2},'/opt/',dakota_filename_in,'.in workdir/',Study_paths{1,1},'/',Study_paths{1,2},'/opt/',dakota_filename_out,'.txt'')');   % unix(''python test_saveFile.py'')
% evalc(shell_command);



for ii = 2:num_studies

    fin  = fopen(strcat('./workdir/',Study_paths{1,1},'/',Study_paths{1,2},'/opt/', dakota_filename_out,'.in'));
    fout = fopen(strcat('./workdir/',Study_paths{ii,1},'/',Study_paths{ii,2},'/opt/', dakota_filename_out,'.in'),'w');
    
    while ~feof(fin)
        f_line = fgetl(fin);
        f_line = strrep(f_line, 'Study0035/0530', strcat(Study_paths{ii,1},'/',Study_paths{ii,2}));
        fprintf(fout,'%s\n',f_line);

    end
    
    fclose(fin);
    fclose(fout);
    
end
%delete(strcat('workdir/',Study_paths{1,1},'/',Study_paths{1,2},'/opt/', dakota_filename_out, '.txt'));