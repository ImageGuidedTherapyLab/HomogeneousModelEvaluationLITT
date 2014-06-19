% This function writes a DAKOTA *.in.* file

function DAKOTA_in_writer ( path, mu_opt, alpha_opt, best_iter, opttype  );

% base_cell = cell (31,1);
% base_cell {1,1} = '20 variables';
% base_cell {8,1} = '1.045000000000000e+03 rho_healthy';
% base_cell {9,1} = '9.000000000000000e-01 anfact_healthy';
% base_cell{10,1} = '6.000000000000000e+00 w_0_healthy';
% base_cell{11,1} = '3.700000000000000e+01 body_temp';
% base_cell{14,1} = '2.100000000000000e+01 probe_init';
% base_cell{15,1} = '8.000000000000000e+03 mu_s_healthy';
% base_cell{16,1} = '3.640000000000000e+03 c_p_healthy';
% base_cell{17,1} = '0.000000000000000e+00 robin_coeff';
% base_cell{18,1} = '3.840000000000000e+03 c_blood_healthy';
% base_cell{22,1} = '2 functions';
% base_cell{23,1} = '1 ASV_1:obj_fn_1';
% base_cell{24,1} = '1 ASV_2:obj_fn_2';
% base_cell{25,1} = '5 derivative_variables';
% base_cell{26,1} = '1 DVV_1:x_displace';
% base_cell{27,1} = '2 DVV_2:mu_eff_healthy';
% base_cell{28,1} = '3 DVV_3:y_displace';
% base_cell{29,1} = '4 DVV_4:z_displace';
% base_cell{30,1} = '5 DVV_5:alpha_healthy';
% base_cell{31,1} = '0 analysis_components';

% Copy the DAKOTA *.in.* file to the main directory
path_base = strcat( 'workdir/',path{1,1},'/',path{1,2},'/opt/');
copy_command = strcat('cp ./', path_base, 'optpp_pds.', opttype, '.in.', num2str(best_iter), ' ./', path_base, 'optpp_pds.LOOCV.in.1');
bash_copy_command = strcat( 'unix(''', copy_command, ''')');
disp(bash_copy_command);
evalc( bash_copy_command);

% Find and replace using Bash's Sed commands
orig_mu = '.*mu_eff_healthy';
new_mu = strcat ( num2str(mu_opt),' mu_eff_healthy');
mu_sed_command = strcat( '''sed "2s/', orig_mu, '/', new_mu, '/" ./',path_base,'optpp_pds.LOOCV.in.1 > ./', path_base,'optpp_pds.LOOCV.in.2' );
bash_mu_command = strcat( 'unix(', mu_sed_command, ''')');
evalc( bash_mu_command);

% sed "3s/.*mu_eff_healthy/289.8162 mu_eff_healthy/" ./optpp_pds.LOOCV.in.1 >optpp_pds.LOOCV.in.2

orig_alpha = '.*alpha_healthy';
new_alpha = strcat ( num2str(alpha_opt),' alpha_healthy');
alpha_sed_command = strcat( '''sed "3s/', orig_alpha, '/', new_alpha, '/" ./', path_base,'optpp_pds.LOOCV.in.2 > ./',path_base,'optpp_pds.LOOCV.in.3' );
bash_alpha_command = strcat( 'unix(', alpha_sed_command, ''')');
evalc( bash_alpha_command);

% Because I'm inept at sed, I'll do this:
mv_command = strcat('mv ./', path_base, 'optpp_pds.LOOCV.in.3 ./', path_base,'optpp_pds.LOOCV.in.1');
bash_mv_command = strcat( 'unix(''', mv_command, ''')');
evalc( bash_mv_command);


%     python_command = strcat( 'unix(''python ./brainsearch.py --param_file ./', param_file, ''')'); % unix(''python test_saveFile.py'')
%     disp(python_command );
%     %evalc(python_command);
% 
% % Find and replace using Bash's Sed commands
% sed_command = strcat( '''sed ''s/', orig_mu, '/', new_mu, '''' );
% 
% 
% fin = fopen(strcat('workdir/',path{1,1},'/',path{1,2},'/opt/', 'optpp_pds.', opttype, '.in.', num2str(best_iter) ) );
% 
% while ~feof(fin)
%     
%     f_line = fgetl(fin);
% 
% 
% 
% fin  = fopen(strcat('workdir/',path{1,1},'/',path{1,2},'/opt/', 'setup.ini'));
% 
% while ~feof(fin)
%         
%     f_line = fgetl(fin);
%     
%     if strcmp (f_line, 'voi = [') == 1
%         base
%     
%     f_line = strrep(f_line, 'Study0035/0530', strcat(Study_paths{ii,1},'/',Study_paths{ii,2}));
%     fprintf(fout,'%s\n',f_line);
%     
% end
end
