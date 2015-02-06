clear
close all
clc

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

data_filename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file

opttype = 'bestfit50' ;
datasummary = dlmread(data_filename,',',1,0);
datasummary(any(isnan(datasummary), 2), 7) = 1;
num_studies = size(datasummary,1);

opt.paths = cell (num_studies,2);
    
Study_paths{1,1} = strcat( 'Study00',num2str(datasummary(1,1)));
Study_paths{1,2} = strcat( '0',num2str(datasummary(1,2)));


% From mu_eff_data, find the matching study's(ies') mu_eff value(s)
total = cell(2,5);
mat_string = ['workdir/',Study_paths{1,1},'/',Study_paths{1,2},'/opt/optpp_pds.',opttype,'.in.1.mat'];
load (mat_string);

n_sources = [1e+0 1e+1 1.5e+1 2.5e+1 5e+1 7.5e+1 1e+2 1e+3 1e+4 1e+5];
%n_sources = [1e+0 1.2e+1 1e+1 1e+2];
total{1,1}{1} = [Study_paths{1,1},'/',Study_paths{1,2}];
total{2,1}{1} = [Study_paths{1,1},'/',Study_paths{1,2}];
total{1,1}{2} = 200;
total{2,1}{2} = 2000;

xx= inputdatavars.voi(2)-inputdatavars.voi(1)+1;
yy= inputdatavars.voi(4)-inputdatavars.voi(3)+1;
l_n_sources = length(n_sources);

for ii = 1:2
    mu = total{ii,1}{2};
    temp_fields = zeros(yy,xx,l_n_sources);
    for jj = 1:l_n_sources
%         % Display run information
%         disp('Start ')
%         disp(strcat (num2str(ii),' of ', num2str(num_studies)))
%         fprintf('iter %s \n', total{ii,1});
        
        % Do global optimization
        [temp_fields(:,:,jj)] = temperature_obj_fxn_GPU_tmap ( inputdatavars, n_sources(jj), mu );
        
    end
    total{ii,2}=temp_fields;
end
clear ii jj
cd ../../../MATLAB/Tests/in_silico/
save ('in_silico_test.mat','total');
cd (path22);

for ii = 1:2
    total{ii,3} = zeros(l_n_sources,3);
    total{ii,4} = zeros( l_n_sources,15);
    total{ii,5} = total{ii,4};
    total{ii,6} = total{ii,4};

    for jj = 1:l_n_sources
        
        tmap_iter = total{ii,2}(:,:,jj);
        tmap_gold = total{ii,2}(:,:,end);
        t_diff = tmap_gold - tmap_iter;

        total{ii,3}(jj,1)= n_sources(jj);
        total{ii,3}(jj,2)= ( norm ( t_diff , 2 ) )^2; % L2 norm
        total{ii,3}(jj,3)= max(max( total{ii,2}(:,:,jj)  ) ) - max( tmap_gold(:));  % 
        
        for kk = 1:15
            model_deg_threshold = tmap_iter >= ( 50 + kk);
            gold_deg_threshold = tmap_gold >= ( 50 + kk);
            n_model = sum( sum( model_deg_threshold ));
            n_gold =  sum( sum( gold_deg_threshold  ));
            intersection = model_deg_threshold + gold_deg_threshold;
            intersection = intersection > 1;
            n_intersection = sum( sum( intersection ));
            total{ii,4}(jj,kk) = n_model - n_intersection; % False positive count
            total{ii,5}(jj,kk) = n_gold - n_intersection; % False negative count
            total{ii,6}(jj,kk) = 2 * n_intersection / ( n_model + n_gold );
        end
        
    end
end

cd ../../../MATLAB/Tests/in_silico/
save ('in_silico_test.mat','total');
cd (path22);


