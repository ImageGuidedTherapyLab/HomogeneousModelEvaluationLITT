clear ii jj
close all
clear

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% cd ../../../MATLAB/Tests/direct_search/
cd ../../../MATLAB/Tests/direct_search/libraries

choice=5;
toss_choice = 1;
metric_choice = 1;  % This is only relevant if choice == 5. 1 is DSC, 2 is L2
if choice == 1   % mu
    
    %load ('GPU_global_mu2.mat');
    load ('GPU_dict_mu.mat');
    
elseif choice == 2  % perf
    
    %load ('GPU_global_perf2.mat');
    load ('GPU_dict_perf.mat');
    
elseif choice == 3   % cond
    
    %load ('GPU_global_cond2.mat');
    load ('GPU_dict_cond.mat');
elseif choice ==5  % random
    
    load ('GPU_dict_perf_mu_rand.mat');
    
    
end

cd(path22);

total(1,:)=[];

if toss_choice == 1
    ix       =find(~cellfun(@isempty,regexp(total(:,1),'0497'))==1); % Absolutely should be excluded
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0378'))==1); % Strongly suggest exclusion
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1); % Strongly suggest exclusion
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1); % Absolutely should be excluded
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0466'))==1); % Very probably suggest exclusion
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0468'))==1); % Very probably suggest exclusion
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0471'))==1); % Strongly suggest exclusion

    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0417'))==1); % Very probably suggest exclusion
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0409'))==1); % Absolutely should be excluded
    
    ix(end+1)=find(~cellfun(@isempty,regexp(total(:,1),'0415'))==1); % Absolutely should be excluded
    total(ix,:) = [];
    
    % ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
    % total(ix,:) = [];
end

aa = cell2mat( total(:,8));

if choice == 1   % mu
    
    index=find( (aa(:,2)>-1)==1);
    
elseif choice == 2  % perf
    
    index=find( (aa(:,2)<26)==1);
    
elseif choice == 3   % cond
    
    index=find( (aa(:,2)>0.4)==1);
    
elseif choice ==5
    index=find( (aa(:,2)>-1)==1);
    
end

if choice ==5
    
    w_lim = [3.25 16];
    %w_lim = [3 16.5];
    w_pix = 250;
    w_array = linspace( w_lim(1),w_lim(2),w_pix);
    mu_lim = [55 400];
    %mu_lim = [50 3000];
    mu_pix = 250;
    mu_array = linspace( mu_lim(1),mu_lim(2),mu_pix);
    %[paraXq, paraYq] = meshgrid (w_lim(1): 0.25 : w_lim(2), mu_lim(1): 50: mu_lim(2) );
    [paraXq, paraYq] = meshgrid ( w_array, mu_array);
    
    %     bb=hist2(total{ii,2}(:,1),total{ii,2}(:,2),200);
    %     figure;imagesc(bb); xlabel('perf');ylabel('mu');
    Xx = zeros(250, 250, length(index));
    Yy = Xx;
    obj_fxn = Xx;
    for jj=1:length(index)
        ii =index(jj);
        
        if metric_choice ==1
            [Xx(:,:,jj), Yy(:,:,jj), obj_fxn(:,:,jj)]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,3}(:,7),paraYq,paraXq);
            %figure; contourf(x,y,z);colorbar; caxis([0 0.9]); title(['DSC ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
            %         elseif metric_choice==2
            %             [x, y, z]=griddata(total{ii,2}(:,1),total{ii,2}(:,2),total{ii,2}(:,4),paraYq,paraXq);
            %             figure; contourf(x,y,z);colorbar; title(['L_2 ' total{ii,1}]); xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]');
        end
        
        
        
    end
    clear ii jj

    opt_mean = mean(obj_fxn,3);
    figure; contourf(Xx(:,:,1),Yy(:,:,1),opt_mean);caxis([0 0.9]); colorbar; title(['Global Mean DSC']);
    xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
    
    opt_median = median(obj_fxn,3);
    figure; contourf(Xx(:,:,1),Yy(:,:,1),opt_median);caxis([0 0.9]); colorbar; title(['Global Median DSC']);
    xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
    
    opt_std = std(obj_fxn,0,3);
    figure; contourf(Xx(:,:,1),Yy(:,:,1),opt_std); caxis([0 0.2]); colorbar; title(['Global St. Dev. of DSC']);
    xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
    
    opt_pass=obj_fxn;
    opt_pass(opt_pass< 0.7)=0;
    opt_pass(opt_pass>=0.7)=1;
    opt_pass = sum(opt_pass,3);
    figure; contourf(Xx(:,:,1),Yy(:,:,1),opt_pass);colorbar; title(['Map of passing datasets']);
    xlabel('mu_{eff}   [ m^{-1} ]'); ylabel('perf [ kg/(m^3 s) ]'); set(findobj('type','axes'),'fontsize',15);
    
    obj_fxn_sz = size(obj_fxn);
    LOOCV_mean = zeros( obj_fxn_sz(1),obj_fxn_sz(2),obj_fxn_sz(3) );
    LOOCV_median = LOOCV_mean;
    LOOCV_std = LOOCV_mean;
    mean_ix1 = zeros( length(index),1);
    mean_ix2 = mean_ix1;
    mean_max = mean_ix1;
    median_ix1 = mean_ix1;
    median_ix2 = mean_ix1;
    median_max = mean_ix1;
    pass_ix1 = mean_ix1;
    pass_ix2 = mean_ix1;
    pass_max = mean_ix1;
    
    for jj=1:length(index)
        ii = index(jj);
        
        obj_fxn_iter = obj_fxn;
        obj_fxn_iter( :,:,ii ) = [];
        LOOCV_mean = mean(obj_fxn_iter,3);
        LOOCV_median = median( obj_fxn_iter,3);
        
        obj_fxn_iter_pass = obj_fxn_iter;
        obj_fxn_iter_pass( obj_fxn_iter_pass < 0.7 ) = 0;
        obj_fxn_iter_pass( obj_fxn_iter_pass >=0.7 ) = 1;
        LOOCV_pass = sum( obj_fxn_iter_pass, 3);
 
        [mean_max(ii), II] = max( LOOCV_mean(:));
        [mean_ix1(ii), mean_ix2(ii)] = ind2sub( size(LOOCV_mean), II );
        
        [median_max(ii), II] = max( LOOCV_median(:));
        [median_ix1(ii), median_ix2(ii)] = ind2sub( size(LOOCV_median), II );
        
        [pass_max(ii), II] = max( LOOCV_pass(:));
        [pass_ix1(ii), pass_ix2(ii)] = ind2sub( size(LOOCV_pass), II );
        
    end
    LOOCV_mean = [ mean_max mean_ix1 mean_ix2];
    LOOCV_median = [ median_max median_ix1 median_ix2];
    LOOCV_pass = [ pass_max pass_ix1 pass_ix2];
             
else
    
    
    
    
    for jj=1:length(index)
        ii = index(jj);
        figure(ii);
        [AX h1 h2] = plotyy(total{ii,2}(:,1),total{ii,2}(:,3),total{ii,2}(:,1),total{ii,3}(:,7));
        legend([h1;h2],'L_2','DSC');
    end
    figure;
    %aa(index,:)=[];
    hist(aa(:,2));
    var=Descriptive_statistics(aa(:,2))
    DSC=Descriptive_statistics_LOOCV(aa(:,1))
    
end



