% This script finds the best mu_eff for the different studies.
function display_var_choice (choice, patient_ix);
%choice = 1; % 1 = mu; 2 = perf; 3 = cond;
tic
close all
% Identify the studies to be examined.
cd /mnt/FUS4/data2/sjfahrenholtz/gitMATLAB/opt_new_database/PlanningValidation
data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );
opttype = 'bestfit51' ;

datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);

opt.paths = cell (num_studies,2);
for ii = 1:num_studies
    
    Study_paths{ii,1} = strcat( 'Study00',num2str(datasummary(ii,1)));
    Study_paths{ii,2} = strcat( '0',num2str(datasummary(ii,2)));
    
end

clear ii
cd ../../../MATLAB/Tests/direct_search

if choice == 1
    
    load ('GPU_global_mu2.mat');
    
elseif choice == 2
    
    load ('GPU_global_perf2.mat');
    
elseif choice == 3
    
    load ('GPU_global_cond2.mat');
    
elseif choice == 4
    cd ../../../MATLAB/Tests/direct_search/libraries
    load ('GPU_dict_perf_mu_global_400');
    
end
total(1,:)=[];
cd (path22);
var_opt = cell2mat( total(:,8));
% indexC = strfind(opt.paths,'Study0035');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom,:)=[];
% datasummary(toss_index_phantom,:)=[];
% num_studies = size(datasummary,1);
%
% indexC = strfind(opt.paths,'0457');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom-num_studies,:)=[];
% datasummary(toss_index_phantom-num_studies,:)=[];
% num_studies = size(datasummary,1);
%
% indexC = strfind(opt.paths,'0476');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom-num_studies,:)=[];
% datasummary(toss_index_phantom-num_studies,:)=[];
% num_studies = size(datasummary,1);
%
% indexC = strfind(opt.paths,'0436');
% toss_index_phantom = find(not(cellfun('isempty',indexC)));
% opt.paths(toss_index_phantom-num_studies,:)=[];
% datasummary(toss_index_phantom-num_studies,:)=[];

% 0491 is suspect: small ROI, ablation fills entire ROI (ie not good)
% 0490 is suspect: small ROI, ablation fills entire ROI
% 0378 is suspect: has very funny shape
% 0385 is really small (but nice and round)
% 0476 is really small (and funny shaped)
% 0435 is really funny shaped
% 0436 small
% 0466 small
% 0468 Really small
% 0471 small and funny shaped
% 0477 small
% 0450 almost small

% From mu_eff_data, find the matching study's(ies') mu_eff value(s)
if choice ==4
    for ii = patient_ix
        path_base = strcat ( 'workdir/',Study_paths{ii,1}, '/', Study_paths{ii,2}, '/opt');
        load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );
        
        [ tmap_model, MRTI_crop, intersection ] = display_obj_fxn_GPU_choice ( inputdatavars, 50, var_opt(ii,2) , var_opt(ii,3), summary.k_cond, choice );
        
        inputdatavars.voi = double( inputdatavars.voi);
        x_diff = abs(inputdatavars.voi(2) - inputdatavars.voi(1));
        x_lim = 0: inputdatavars.spacing(1): x_diff*inputdatavars.spacing(1);
        y_diff = abs(inputdatavars.voi(4) -inputdatavars.voi(3));
        y_lim = 0: inputdatavars.spacing(2): y_diff*inputdatavars.spacing(2);
        
        x_lim(end) % FOV limits
        y_lim(end)
        
        figure; imagesc(tmap_model, [30 100]); colorbar; set(findobj('type','axes'),'fontsize',14);
        figure; imagesc(MRTI_crop, [30 100]); colorbar; set(findobj('type','axes'),'fontsize',14);
        figure; imagesc(intersection(:,:,7)); title( Study_paths{ii,2} ); colorbar; set(findobj('type','axes'),'fontsize',14);
        
        aa = cell2mat( total(:,8));
        index=find( (aa(:,2)>-1)==1);
        w_array = unique( summary.w_perf);
        mu_array = unique( summary.mu);
        [paraXq, paraYq] = meshgrid ( w_array, mu_array);
        
        %     bb=hist2(total{ii,2}(:,1),total{ii,2}(:,2),200);
        %     figure;imagesc(bb); xlabel('perf');ylabel('mu');
        grid_sz = size (paraXq);
        Xx = zeros(grid_sz(1), grid_sz(2), length(index));
        Yy = Xx;
        obj_fxn = Xx;
        
        [paraXq, paraYq] = meshgrid ( w_array, mu_array);
        
        %     bb=hist2(total{ii,2}(:,1),total{ii,2}(:,2),200);
        %     figure;imagesc(bb); xlabel('perf');ylabel('mu');
        %         for jj=1:length(index)
        %             kk =index(jj);
        dice = squeeze(total{ii,3});
        
        %[Xx(:,:,jj), Yy(:,:,jj), obj_fxn(:,:,jj)]=griddata(total{kk,2}(:,1),total{kk,2}(:,2),dice,paraYq,paraXq);
        [Xx, Yy, obj_fxn]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),dice,paraYq,paraXq);
        
        %end
        clear kk
        
        figure; contourf(Xx,Yy,obj_fxn);caxis([0 0.9]); colorbar; title( Study_paths{ii,2} );
        xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
        var_opt(ii,:)
        %[distXq, distYq] = meshgrid (x_lim, y_lim);
%         figure; contourf(distXq,distYq,intersection(:,:,7), [0 1 2 3]); colorbar; title( Study_paths{ii,2} );
%         xlabel('Distance   [ m ]'); ylabel('Distance [ m ]'); set(findobj('type','axes'),'fontsize',15);
        keyboard
        close all
    end
    
else
    
    for ii = patient_ix
        path_base = strcat ( 'workdir/',Study_paths{ii,1}, '/', Study_paths{ii,2}, '/opt');
        load( strcat ( path_base, '/optpp_pds.', opttype, '.in.1.mat') );
        %inputdatavars.spacing
        if choice == 1  % mu
            
            [ tmap_model, MRTI_crop, intersection ] = display_obj_fxn_GPU_choice ( inputdatavars, 25, var_opt(ii,2) , summary.w_perf, summary.k_cond, choice );
            
        elseif choice == 2 % perf
            
            [ tmap_model, MRTI_crop, intersection ] = display_obj_fxn_GPU_choice ( inputdatavars, 25, summary.mu_eff , var_opt(ii,2), summary.k_cond, choice );
            
        elseif choice == 3 % cond
            
            [ tmap_model, MRTI_crop, intersection ] = display_obj_fxn_GPU_choice ( inputdatavars, 25, summary.mu_eff , summary.w_perf, var_opt(ii,2), choice );
            
        end
        inputdatavars.voi = double( inputdatavars.voi);
        x_diff = abs(inputdatavars.voi(2) - inputdatavars.voi(1));
        x_lim = 0: inputdatavars.spacing(1): x_diff*inputdatavars.spacing(1);
        y_diff = abs(inputdatavars.voi(4) -inputdatavars.voi(3));
        y_lim = 0: inputdatavars.spacing(2): y_diff*inputdatavars.spacing(2);
        
        x_lim(end)
        y_lim(end)
        
        figure; imagesc(tmap_model, [30 100]); colorbar; set(findobj('type','axes'),'fontsize',14);
        figure; imagesc(MRTI_crop, [30 100]); colorbar; set(findobj('type','axes'),'fontsize',14);
        figure; imagesc(intersection(:,:,7)); title( Study_paths{ii,2} ); colorbar; set(findobj('type','axes'),'fontsize',14);
        figure; [AX h1 h2] = plotyy(total{ii,2}(:,1),total{ii,2}(:,3),total{ii,2}(:,1),total{ii,3}(:,7));
        set(h1, 'LineWidth',5); set(h2,'LineWidth',5); set(findobj('type','axes'),'fontsize',16);
        legend([h1;h2],'L_2','DSC'); title( Study_paths{ii,2});
        var_opt(ii,:)
        keyboard
        close all
        
        
    end
    clear ii
end
end
