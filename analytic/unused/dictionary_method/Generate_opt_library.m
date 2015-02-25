function [ all_opt_fig, no_pwr_fig, sim_dim, summary ] = Generate_opt_library ( max_phys_sz, choice, Study_paths, opttype );

% Make the LOOCV iteration system

path_base = strcat ( 'workdir/',Study_paths{1,1}, '/', Study_paths{1,2}, '/opt');
load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );
mu = str2num( inputdatavars.cv.mu_eff_healthy );
kk = inputdatavars.cv.k_0;
ww = inputdatavars.cv.w_0;


if choice ==1   % Mu
    mu_eff(1) = 0.008;
            % mu_eff(2:101) = linspace(101,200,100);
    %mu_eff=[200 300 400 500 600 700];
    %mu_eff(2:1001) = linspace(1,1000,1000);
    %mu_eff=180;
    %mu_eff(2:301) = linspace(1,300,300);
    mu_eff(2:5001) = linspace(1,5000,5000);
    %mu_eff(1) = 0.008;
    %mu_eff(2:101) = linspace(101,200,100);
    %mu_eff = linspace(100,500,401);
    %mu_eff = [200 500];
    %k_cond = 0.527;
    %w_perf = 6;
    k_cond = kk;
    w_perf = ww;

elseif choice ==2   % perf
    %w_perf = 0.01;
    %w_perf(2:35) = linspace ( 0.5, 17, 34);
    %w_perf(2:34) = linspace ( 0.5, 16.5,33);
    %w_perf(2:133) = linspace ( 0.125, 16.5,132);
    w_perf = linspace( 0.5, 16.5, 129);
    
    % mu_eff = 180;
    % k_cond = 0.527;
    mu_eff = 180;
    k_cond = kk;
    
elseif choice ==3    % cond
    %k_cond = linspace ( 0.01,2,200);
    k_cond = linspace ( 0.52, 3, 249);
    %k_cond = linspace ( 0.52, 3, 2);
    k_cond(end+1) = 0.527;
    k_cond = sort(k_cond);
    %k_cond = [0.527 0.55];
    % mu_eff = 180;
    % w_perf = 6;
    mu_eff = 180;
    w_perf = ww;
    
elseif choice ==4 % perf and mu
    %w_perf = linspace( 0.5, 16.5, 129);
    w_perf = [6 9 12];
    mu_eff =[200 300 400 500 600 700];
    k_cond = kk;
   
elseif choice ==5 % perf and mu randomly sampled
    
    w_perf = 0.5 : 0.125 : 15;
    mu_eff = 50:4500;
    w_Num_runs = 50;
    m_Num_runs = 1000;
    Num_runs = w_Num_runs * m_Num_runs;
    param_list = zeros ( Num_runs,3);
    w_list = sort( datasample( w_perf,w_Num_runs,'Replace',true) );
    mu_list = sort( datasample( mu_eff, m_Num_runs,'Replace',true) );
    param_list(:,1) = 1:Num_runs;
    %for ii = 1:w_Num_runs
        
        
    
    
    w_perf = linspace( 0.5, 16.5, 129);
    w_perf = [6 9 12];
    mu_eff =[200 300 400 500 600 700];
    k_cond = kk;
    
end
    summary.k_cond = k_cond;
    summary.w_perf = w_perf;
    summary.mu = mu_eff;
    summary.choice = choice;
    [all_opt_fig, no_pwr_fig,sim_dim] = temperature_GPU_1W ( inputdatavars, 50, max_phys_sz, mu_eff, w_perf, k_cond, choice );
end