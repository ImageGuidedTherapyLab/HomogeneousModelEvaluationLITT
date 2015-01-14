% This script finds all of the optimized outputs

clear all
cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds/
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );
Patient_Paths = importdata( 'patient_paths.txt' );
obj_fxn = zeros(28,2);
power_log = zeros(28,1);
mu_eff = zeros(28,1);
keep_index = 1;
keep_fxn = zeros(2,1);
%keep_patient_path (1:2,1) = 'a';


for patient_index = 1:28
    one_patient_path = char ( Patient_Paths ( patient_index ) );
    patient_opt_path = strcat( path22, one_patient_path);
    cd (patient_opt_path);
    tabular_dat = load ( 'dakota_q_newton_heating.in.tabular.dat' );
    null_check = isempty (tabular_dat);
    if null_check == 1;
        tabular_dat = 0;
    end
    [obj_fxn( patient_index , 1 ), obj_fxn( patient_index , 2 )] = min( tabular_dat ( :,3) );
    
    input_param = 'optpp_pds.in.';
    index = num2str( obj_fxn( patient_index , 2 ) );
    input_filename = strcat( input_param, index);
    input_filename = strcat( input_filename, '.mat' );
    
    load(input_filename);
    mu_eff ( patient_index) = str2num(mu_eff_healthy);
    
    power_one_patient = load ('time_then_power.csv');
    power_log (patient_index) = max( power_one_patient (:,2));
    clear power_one_patient
    
    if obj_fxn( patient_index,1) < 3e+5
        keep_fxn ( keep_index) = obj_fxn( patient_index, 1);
        keep_mu_eff ( keep_index) = mu_eff( patient_index);
        keep_patient_path ( keep_index) = Patient_Paths(patient_index);
        keep_index = keep_index + 1;
    end
    
    
    
    clear obj_fxn_list
    

end

clear anfact_healthy input_filename input_param k_0_healthy mu_a_healthy mu_s_healthy one_patient_path path22 patient_index patient_opt_path probe_init
clear robin_coeff w_0_healthy x_displace x_rotate y_displace y_rotate z_displace z_rotate

% 
% ii = 1;
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0002/000/opt
% load ( 'obj_fxn_list.txt' );
% obj_fxn( ii ) = obj_fxn_list( ii );
% clear obj_fxn_list
% ii = ii + 1;
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0002/001/opt
% load ( 'obj_fxn_list.txt' );
% obj_fxn( ii ) = obj_fxn_list( ii );
% clear obj_fxn_list
% ii = ii + 1;
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0002/010/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0002/017/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0002/020/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0002/021/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0002/022/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0003/002/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0003/013/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0003/018/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0003/025/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0004/003/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0004/006/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0005/004/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0005/008/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0005/012/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0006/005/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0006/007/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0006/009/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0006/019/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0006/023/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0006/024/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0007/011/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0007/015/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0007/026/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0008/014/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0008/016/opt
% 
% cd /FUS4/data2/sjfahrenholtz/gitMATLAB/optpp_pds
% cd workdir/Patient0008/027/opt
