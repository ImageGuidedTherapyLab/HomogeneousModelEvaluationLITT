clear ii jj
close all
clear

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

cd ../../../MATLAB/Tests/direct_search

choice=3;
if choice == 1   % mu
    
    load ('GPU_global_mu.mat');
    
elseif choice == 2  % perf
    
    load ('GPU_global_perf.mat');
    
elseif choice == 3   % cond
    
    load ('GPU_global_cond.mat');
    
end

cd(path22);

total(1,:)=[];

aa = cell2mat( total(:,8));

if choice == 1   % mu
    
    index=find( (aa(:,2)>700)==1);
    
elseif choice == 2  % perf
    
    index=find( (aa(:,2)>6)==1);
    
elseif choice == 3   % cond
    
    index=find( (aa(:,2)>0.75)==1);
    
end

for jj=1:length(index)
    ii = index(jj);
    figure(ii);
    [AX h1 h2] = plotyy(total{ii,2}(:,1),total{ii,2}(:,3),total{ii,2}(:,1),total{ii,3}(:,7));
    legend([h1;h2],'L_2','DSC');
end