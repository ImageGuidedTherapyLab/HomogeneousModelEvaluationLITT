clear ii jj
close all
clear

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

% cd ../../../MATLAB/Tests/direct_search/
cd ../../../MATLAB/Tests/direct_search/libraries

choice=2;
if choice == 1   % mu
    
    %load ('GPU_global_mu2.mat');
    load ('GPU_dict_mu.mat');
    
elseif choice == 2  % perf
    
    %load ('GPU_global_perf2.mat');
    load ('GPU_dict_perf.mat');
    
elseif choice == 3   % cond
    
    %load ('GPU_global_cond2.mat');
    load ('GPU_dict_cond.mat');
    
end

cd(path22);

total(1,:)=[];

% List of excluded datasets   
total(1,:) = []; % Drop the labels
ix=find(~cellfun(@isempty,regexp(total(:,1),'0497'))==1); % Absolutely should be excluded
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0378'))==1); % Strongly suggest exclusion
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0476'))==1); % Strongly suggest exclusion
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0436'))==1); % Absolutely should be excluded
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0466'))==1); % Very probably suggest exclusion
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0468'))==1); % Very probably suggest exclusion
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0471'))==1); % Strongly suggest exclusion
total(ix,:) = [];
% 
% ix=find(~cellfun(@isempty,regexp(total(:,1),'0417'))==1); % Very probably suggest exclusion
% total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0409'))==1); % Absolutely should be excluded
total(ix,:) = [];

ix=find(~cellfun(@isempty,regexp(total(:,1),'0415'))==1); % Absolutely should be excluded
total(ix,:) = [];


% ix=find(~cellfun(@isempty,regexp(total(:,1),'0457'))==1);
% total(ix,:) = [];


aa = cell2mat( total(:,8));

if choice == 1   % mu
    
    index=find( (aa(:,2)>-1)==1);
    
elseif choice == 2  % perf
    
    index=find( (aa(:,2)<26)==1);
    
elseif choice == 3   % cond
    
    index=find( (aa(:,2)>0.4)==1);
    
end

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